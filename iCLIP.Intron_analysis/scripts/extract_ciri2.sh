set -e

read circRNA_tab circRNA_reference <<< "$@"
head -n 1 $circRNA_reference | grep "^circRNA_ID"
grep -f <(grep -v circRNA_ID $circRNA_tab) $circRNA_reference
