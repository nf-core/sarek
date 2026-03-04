# Plot Generation in Rust

Guidance for generating publication-quality plots that match upstream tool output. Most
bioinformatics tools produce plots (QC summaries, coverage profiles, density scatters),
and matching them exactly is often the most tedious phase of a rewrite.

## General Principles

- **Budget significant time.** Plot generation and visual matching can consume 1–3 sessions
  per tool. It's iterative and cannot be fully automated.
- **Match upstream default figure dimensions.** R defaults to roughly square plots (480×480
  or similar). Using 640×480 changes aspect ratios of all elements. Check the upstream source
  for explicit `width`/`height` settings.
- **Always generate both SVG and PNG.** SVG is more portable for web, docs, and MultiQC
  integration. PNG is needed for pixel-level comparison and reports that don't support SVG.
- **Match color palettes and parameters exactly.** Copy color hex values, density estimation
  bandwidths, bin counts, and axis ranges directly from the upstream source. Never assume
  from memory — always verify by running the upstream tool and inspecting actual output.
  (e.g., R's `densCols` default ramp starts at **black**, not cyan.)
- **Establish a fast visual comparison pipeline early.** Automate: build → regenerate plots →
  side-by-side A/B comparison (e.g., open reference and generated images together). Track
  progress with quality scores (1–10 scale) per plot.
- **Don't hand-roll complex SVG charts.** Use a proper charting library (`plotters` for Rust,
  `matplotlib` for Python helper scripts). Hand-written SVG for things like benchmark charts
  is fragile and will have label positioning issues.

## Plotters Crate: Known Issues and Workarounds

The `plotters` crate (v0.3) is the most mature Rust charting library but has significant
gaps for scientific plotting. Expect to work around these:

### Rotated Axis Labels Are Broken

Disable automatic labels, then draw manually:

```rust
// Disable automatic x-axis labels
chart.configure_mesh()
    .x_label_formatter(&|_| String::new())
    .draw()?;

// Draw rotated labels manually
for (value, label) in tick_values.iter() {
    let (px, py) = chart.plotting_area().map_coordinate(&(*value, y_min));
    root.draw_text(
        label,
        &text_style.transform(FontTransform::Rotate270),
        (px - 3, py + 5),  // Adjust offset to align with tick marks
    )?;
}
```

### `y_labels(N)` Doesn't Give N Labels

Plotters picks "nice" intervals, so `y_labels(5)` on a 0–100 range may produce just
0/50/100 instead of 0/25/50/75/100. Fix: request many more labels and filter in the
formatter:

```rust
chart.configure_mesh()
    .y_labels(21)  // Request many labels (step=5 on 0–100)
    .y_label_formatter(&|v| {
        if (*v as i32) % 25 == 0 { format!("{}", *v as i32) }
        else { String::new() }
    })
    .draw()?;
```

### Thick Lines Look Blocky (No Anti-Aliasing)

`LineSeries` with `stroke_width > 1` produces jagged, aliased lines with visible step
artifacts. This is especially noticeable on dashed curves (fit lines, model overlays).

**Workaround: Draw curves as filled circles along the path.**

Instead of a thick `LineSeries`, plot small filled `Circle` primitives placed at regular
intervals along the curve using arc-length parameterization:

1. Convert data points to pixel coordinates via `chart.plotting_area().map_coordinate()`
2. Walk along the pixel polyline by arc-length, accumulating distance
3. At each step, interpolate the data coordinates at fractional positions within segments
4. Alternate between "dash" segments (draw circles) and "gap" segments (skip)
5. Draw all circles via `chart.draw_series()` so they clip to the plot area

This produces dramatically smoother curves than `LineSeries`, especially for dashed lines.
The `radius` parameter controls line thickness, and `dash_px`/`gap_px` control dash pattern.

### No Built-In Pie Chart

Render manually using `Polygon::new` for each slice:

```rust
// For each slice, compute start/end angles and create a polygon
let segments = 50; // Arc resolution
for (start_angle, end_angle, color) in slices {
    let mut points = vec![(cx, cy)]; // Center point
    for i in 0..=segments {
        let angle = start_angle + (end_angle - start_angle) * (i as f64 / segments as f64);
        let px = cx + (radius as f64 * angle.cos()) as i32;
        let py = cy - (radius as f64 * angle.sin()) as i32;
        points.push((px, py));
    }
    root.draw(&Polygon::new(points, color.filled()))?;
}
```

### No Built-In Boxplot

Render from primitives:

- **Box**: `Rectangle::new` with fill color + border stroke
- **Median line**: `LineSeries` (thicker `stroke_width`)
- **Whiskers**: `LineSeries` (thin lines from Q1→lower fence, Q3→upper fence)
- **Caps**: Short horizontal `LineSeries` at fence endpoints
- **Outliers**: `Circle::new` (open circles — use stroke only, not filled)

### Coordinate Mapping

To position custom-drawn elements relative to data coordinates:

```rust
let (px, py) = chart.plotting_area().map_coordinate(&(x_data, y_data));
// (px, py) is a BackendCoord (i32, i32) — pixel position
```

This is essential for drawing manual annotations, custom legends, or the circle-based
curves described above.

## Resolution-Aware Scaling

Create a scaling helper that keeps PNG and SVG rendering consistent:

```rust
const SCALE: f64 = 4.0; // PNG resolution multiplier

fn render_plot<DB: DrawingBackend>(root: DrawingArea<DB, Shift>, data: &Data, pxs: f64) {
    let ps = |v: f64| (v * pxs) as u32;

    // Use ps() for dimensions that should scale:
    let font_size = ps(12.0);     // Font size scales linearly
    let margin = ps(10.0);        // Margins scale linearly
    let stroke = ps(1.0);         // Line width scales linearly

    // Some dimensions should NOT scale linearly:
    let point_radius = (pxs * 1.2).round().max(1.0) as u32;  // Sub-linear for scatter
    let min_line_width = 1;  // Never thinner than 1px
    // ...
}

// Public API generates both formats from the same rendering code:
pub fn plot_chart(path: &Path, data: &Data) -> Result<()> {
    let (w, h) = (480, 480);  // Base dimensions matching R defaults

    // PNG at high resolution
    let png_root = BitMapBackend::new(&png_path, (w * SCALE as u32, h * SCALE as u32)).into_drawing_area();
    render_plot(png_root, data, SCALE)?;

    // SVG at base resolution
    let svg_root = SVGBackend::new(&svg_path, (w, h)).into_drawing_area();
    render_plot(svg_root, data, 1.0)?;
    Ok(())
}
```

## Scatter Plot Optimization

### Pixel Deduplication

When plotting thousands of data points, many map to the same pixel. Deduplicate to reduce
SVG file size and speed rendering:

```rust
let mut pixel_map: HashMap<(i32, i32), (f64, f64, f64)> = HashMap::new(); // pixel -> (x, y, density)
for &(x, y, density) in &points {
    let px = chart.plotting_area().map_coordinate(&(x, y));
    pixel_map.entry(px)
        .and_modify(|e| { if density > e.2 { *e = (x, y, density); } })
        .or_insert((x, y, density));
}
```

### Draw Order

Sort points ascending by importance (e.g., density) so that high-value points render on
top and aren't obscured by low-value points.

## Density Estimation

If the upstream tool uses R's `densCols()` / `bkde2D` (2D kernel density estimation),
you need to match the algorithm:

- **Bandwidth selection**: R uses `KernSmooth::dpik()` which implements a plug-in method.
  Silverman's rule (`0.9 * min(sd, IQR/1.34) * n^(-1/5)`) is a reasonable approximation.
- **Separable Gaussian smoothing**: Apply 1D Gaussian smoothing in X, then Y (or vice
  versa). This is mathematically equivalent to 2D Gaussian but much faster.
- **Kernel truncation**: R's `bkde2D` uses `tau = 3.4` standard deviations by default.
  Points beyond this distance contribute zero.
- **Grid resolution**: R defaults to a grid (e.g., 500 bins per axis). Match this.
- **Anisotropic bandwidth**: X and Y bandwidths are typically different. R might use
  ~19 bins in X and ~1 bin in Y depending on data spread. Don't assume isotropic.

## Matching R Plot Parameters

When the upstream tool uses R's base `plot()` or `ggplot2`:

| R Parameter | What It Controls | How to Match |
|-------------|-----------------|--------------|
| `cex` | Point size multiplier | Scale point radius: `cex=0.25` → 1px at base resolution |
| `pch` | Point shape (20=filled circle, 1=open circle) | `Circle::new` filled or stroke-only |
| `lwd` | Line width | `stroke_width(ps(lwd))` |
| `lty` | Line type (1=solid, 2=dashed, 3=dotted) | `LineSeries` or dashed-circle technique |
| `col` | Color (name or hex) | Convert R color names to hex; verify with `col2rgb("name")` |
| `xlim`/`ylim` | Axis range | Set via `.x_label_area_size()` and chart range |
| `main` | Title | Draw manually with `Text::new` for precise positioning |
| `las=2` | Rotated labels | Use manual `FontTransform::Rotate270` workaround |

R's `lty=3` means **dotted** (not dashed) — this is a common mistake.

## Common Pitfalls

1. **Color palette wrong from memory.** R's `densCols()` default starts at **black**, not
   cyan. Always run the upstream tool and extract actual hex values.
2. **Rectangular instead of square.** If you use 640×480 when R uses ~480×480, everything
   looks stretched. Match dimensions first.
3. **Font rendering differs.** Plotters uses system fonts via fontconfig. The exact font
   available may differ between Linux/macOS, causing slight layout differences. This is
   usually acceptable.
4. **Point size scaling.** For scatter plots with thousands of points, point radius should
   often be fixed at 1px regardless of scale factor. Sub-linear scaling
   (`(pxs * 1.2).round().max(1.0)`) keeps points visible without becoming blobs.
5. **Legend positioning.** Plotters' built-in legend may not match R's placement. For
   complex legends, draw manually with `PathElement` (line style) + `Circle` (point marker)
   + `Text::new`.
6. **Axis padding.** R auto-pads axes by ~4% beyond data range. If your axis range is
   exactly `[min, max]`, labels at the boundaries may clip. Add padding:
   `let pad = (max - min) * 0.04;`
