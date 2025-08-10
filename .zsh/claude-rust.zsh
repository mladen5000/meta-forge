export CLAUDE_CMD=claude

alias cc='scripts/cc.sh'
alias ccdiff='scripts/cc_diff.sh'
alias ccfix='scripts/cc_fix.sh'

alias cg='cargo'
alias cgt='cargo test'
alias cgc='cargo check'
alias cgl='cargo clippy --all-targets --all-features -- #-D warnings'
alias cgfmt='cargo fmt --all'

alias gsa='git status -sb'
alias gap='git apply --reject --whitespace=fix'
alias gca='git add -A && git commit -m "apply surgical patch"'
