# CRAN Submission Checklist

## Before Submission

- [ ] Run R CMD check --as-cran
- [ ] Check all URLs are valid
- [ ] Ensure version number is updated
- [ ] Update NEWS.md
- [ ] Check spelling with spelling::spell_check_package()
- [ ] Check with rhub::check_for_cran()

## Submission

1. Upload to CRAN: https://cran.r-project.org/submit.html
2. Wait for email confirmation
3. Address any feedback

## Files Included

- OmniSciKit_0.1.0.tar.gz (main package)
- cran-comments.md (check results)
- NEWS.md (version history)
