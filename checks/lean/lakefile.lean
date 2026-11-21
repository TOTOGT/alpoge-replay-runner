import Lake
open Lake DSL

-- Minimal project: pins Mathlib to the same tag as lean-toolchain and lake-manifest.json.
package s3zeta

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.32.0"

lean_lib S3Zeta where
  srcDir := ".."
  roots := #[`S3Zeta]

lean_lib S3Audit where
  roots := #[`S3Audit]
