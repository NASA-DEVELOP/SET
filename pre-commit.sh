STASH_NAME="pre-commit-$(date +%s)"
git stash save -q --keep-index $STASH_NAME

python -m unittest test

STASHES=$(git stash list)
if [[ $STASHES == "$STASH_NAME" ]]; then
	git stash pop -q
fi
