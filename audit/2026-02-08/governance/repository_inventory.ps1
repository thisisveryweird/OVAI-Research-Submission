# Read-only repository evidence; stdout only. No imports of repository code.
param([string]$Root = (Get-Location).Path)
$ErrorActionPreference = 'Stop'
$r = (Resolve-Path -LiteralPath $Root).Path.TrimEnd('\','/')
Get-ChildItem -LiteralPath $r -Force -Recurse -File | Where-Object { $_.FullName -notlike ($r+'\audit\*') } | Sort-Object FullName | ForEach-Object {
 [pscustomobject]@{path=$_.FullName.Substring($r.Length+1).Replace('\','/'); bytes=$_.Length; sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $_.FullName).Hash.ToLower()}
} | ConvertTo-Json -Depth 3
# .git included; all audit outputs excluded. Hashes establish bytes, not validity.
