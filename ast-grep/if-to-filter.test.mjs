import assert from "node:assert/strict";
import { spawnSync } from "node:child_process";
import { mkdirSync, mkdtempSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join, relative, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import test from "node:test";

test("if-to-filter honours source paths, suppression and warning/error severity", () => {
    const root = mkdtempSync(join(tmpdir(), "sarek-structural-check-"));
    const source = "workflow { if (params.align) { ALIGN(reads) } }";
    const included = ["main.nf", "workflows/align.nf", "subworkflows/local/align/main.nf"];
    const excluded = [
        "modules/local/align/main.nf",
        "modules/nf-core/align/main.nf",
        "subworkflows/nf-core/align/main.nf",
        "tests/align.nf",
        "work/align.nf",
        ".nf-test/align.nf",
        "ast-grep/fixture.nf",
        "nextflow.config",
    ];
    try {
        for (const path of [...included, ...excluded]) {
            const target = join(root, path);
            mkdirSync(dirname(target), { recursive: true });
            writeFileSync(target, source);
        }
        writeFileSync(
            join(root, "workflows/suppressed.nf"),
            "workflow {\n// ast-grep-ignore: if-to-filter\nif (params.align) { ALIGN(reads) }\n}",
        );
        for (const scanPath of [".", root]) {
            for (const strict of [false, true]) {
                const result = spawnSync(
                    fileURLToPath(new URL("node_modules/.bin/ast-grep", import.meta.url)),
                    [
                        "scan",
                        "-c",
                        fileURLToPath(new URL("sgconfig.yml", import.meta.url)),
                        "--json=compact",
                        "--no-ignore=hidden",
                        ...(strict ? ["--error=if-to-filter"] : []),
                        scanPath,
                    ],
                    { cwd: root, encoding: "utf8" },
                );
                assert.equal(result.status, strict ? 1 : 0, result.stderr);
                const matches = JSON.parse(result.stdout);
                assert.deepEqual(
                    matches.map((match) => relative(root, resolve(root, match.file))).sort(),
                    [...included].sort(),
                );
                for (const match of matches) {
                    assert.equal(match.ruleId, "if-to-filter");
                    assert.equal(match.severity, strict ? "error" : "warning");
                    assert.equal(match.replacement, undefined);
                }
            }
        }
    } finally {
        rmSync(root, { recursive: true, force: true });
    }
});
