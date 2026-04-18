## Tests for .bw_path_for_rtracklayer()
##
## The helper rewrites Windows drive-letter paths to file:///C:/... URIs
## so the kent-library protocol parser (udcProtNew) does not misread the
## drive letter as a URL scheme. On non-Windows platforms the path is
## passed through unchanged.

test_that(".bw_path_for_rtracklayer is a pass-through on non-Windows", {
  skip_on_os("windows")
  p <- "/tmp/example.bw"
  expect_identical(epiRomics:::.bw_path_for_rtracklayer(p), p)
})

test_that(".bw_path_for_rtracklayer rewrites drive-letter paths on Windows", {
  skip_on_os(c("mac", "linux", "solaris"))
  ## On Windows, normalizePath will coerce slashes; ensure the output
  ## begins with "file:///" and preserves the drive letter.
  p_in <- "C:/Users/test/example.bw"
  out  <- epiRomics:::.bw_path_for_rtracklayer(p_in)
  expect_match(out, "^file:///[A-Za-z]:/")
})

test_that(".bw_path_for_rtracklayer preserves spaces in drive-letter paths", {
  skip_on_os(c("mac", "linux", "solaris"))
  ## Paths with spaces (common under C:/Users/<name with space>/) must
  ## round-trip into the file:/// URI unchanged -- spaces are legal in
  ## file URIs and rtracklayer accepts them.
  p_in <- "C:/Users/Some User/example.bw"
  out  <- epiRomics:::.bw_path_for_rtracklayer(p_in)
  expect_match(out, "^file:///[A-Za-z]:/Users/Some User/example\\.bw$")
})

test_that(".bw_path_for_rtracklayer returns URLs identically on Windows", {
  skip_on_os(c("mac", "linux", "solaris"))
  ## Pre-qualified URIs must be returned byte-for-byte unchanged. The
  ## helper short-circuits any scheme://-prefixed input before
  ## normalizePath runs, so http://, https://, file:///, ftp://, s3://
  ## all pass through.
  urls <- c(
    "http://example.org/x.bw",
    "https://example.org/x.bw",
    "file:///C:/x.bw",
    "ftp://example.org/x.bw",
    "s3://bucket/x.bw")
  for (u in urls) {
    out <- epiRomics:::.bw_path_for_rtracklayer(u)
    expect_identical(out, u)
  }
})

test_that(".bw_path_for_rtracklayer leaves UNC paths unchanged on Windows", {
  skip_on_os(c("mac", "linux", "solaris"))
  ## UNC paths (\\server\share\file.bw) do not match the drive-letter
  ## regex after normalizePath (which yields //server/share/...), so the
  ## helper returns them unmodified. rtracklayer / kent-library handle
  ## UNC paths directly on Windows.
  p_in <- "//server/share/example.bw"
  out  <- epiRomics:::.bw_path_for_rtracklayer(p_in)
  expect_false(base::startsWith(out, "file:///"))
})
