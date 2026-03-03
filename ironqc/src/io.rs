//! Filesystem helpers for deterministic placeholder output creation.

use anyhow::{Context, Result};
use std::fs::{self, File};
use std::path::Path;

/// Create an empty file at `path`, including missing parent directories.
pub fn touch(path: impl AsRef<Path>) -> Result<()> {
    let path = path.as_ref();

    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create parent directory for {}", path.display()))?;
    }

    File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    Ok(())
}

/// Ensure a directory exists.
pub fn ensure_dir(path: impl AsRef<Path>) -> Result<()> {
    let path = path.as_ref();
    fs::create_dir_all(path)
        .with_context(|| format!("failed to create directory {}", path.display()))?;
    Ok(())
}
