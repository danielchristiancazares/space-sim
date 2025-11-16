# Development workflow
- Create a topic branch from main - e.g. feat/interactive-prompt.
- Keep your changes focused. Multiple unrelated fixes should be opened as separate PRs.
- Ensure your change is free of lint warnings and test failures.
# Writing code changes
- Add or update tests. Every new feature or bug-fix should come with test coverage that fails before your change and passes afterwards. 100% coverage is not required, but aim for meaningful assertions.
- Document behaviour. If your change affects behaviour, update the README, documentation, or comments.
- Keep commits atomic. Each commit should compile and the tests should pass. This makes reviews and potential rollbacks easier.
- Run all checks locally (cargo test && cargo clippy --tests && cargo fmt). CI failures that could have been caught locally slow down the process.
- Make sure your branch is up-to-date with main and that you have resolved merge conflicts.
# Making commits
- Ensure CHANGELOG.md is updated following keepachangelog.com
Types of changes
```
Added for new features.
Changed for changes in existing functionality.
Deprecated for soon-to-be removed features.
Removed for now removed features.
Fixed for any bug fixes.
Security in case of vulnerabilities.
- The commit message should be structured as follows:
```
<type>[optional scope]: <description>
```
Where type should be one of the following:
```
build: Changes that affect the build system or external dependencies (example scopes: gulp, broccoli, npm)
ci: Changes to our CI configuration files and scripts (example scopes: Travis, Circle, BrowserStack, SauceLabs)
docs: Documentation only changes
feat: A new feature
fix: A bug fix
perf: A code change that improves performance
refactor: A code change that neither fixes a bug nor adds a feature
style: Changes that do not affect the meaning of the code (white-space, formatting, missing semi-colons, etc)
test: Adding missing tests or correcting existing tests
```
