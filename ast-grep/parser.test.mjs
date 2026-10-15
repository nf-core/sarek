import assert from "node:assert/strict";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import test from "node:test";

test("Nextflow parser matches channel calls and captures metavariables", () => {
    const result = spawnSync(
        fileURLToPath(new URL("node_modules/.bin/ast-grep", import.meta.url)),
        ["run", "-c", "sgconfig.yml", "-l", "nextflow", "-p", "channel.of($VALUE)", "--stdin", "--json=compact"],
        {
            cwd: fileURLToPath(new URL(".", import.meta.url)),
            input: "workflow { channel.of(42).view() }",
            encoding: "utf8",
        },
    );
    assert.equal(result.status, 0, result.stderr);
    const matches = JSON.parse(result.stdout);
    assert.equal(matches.length, 1);
    assert.equal(matches[0].text, "channel.of(42)");
    assert.equal(matches[0].metaVariables.single.VALUE.text, "42");
});
