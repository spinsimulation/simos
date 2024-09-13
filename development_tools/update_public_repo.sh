#!/bin/bash

# Ensure the script exits on any error
set -e

# Ensure the script is running from the root of the repository
if [ ! -d ".git" ]; then
    echo "This script must be run from the root of the repository."
    exit 1
fi

# Ensure that remote repositories are set up
if [ -z "$(git remote | grep github)" ]; then
    echo "The remote repository 'github' is not set up."
    exit 1
fi

# Define the branch names
PRIVATE_BRANCH="ld2"
PUBLIC_BRANCH="public"

# Fetch the latest changes from the public repository
echo "Fetching the latest changes from the public repository..."
echo "[git fetch github]"
git fetch github

# Checkout the private branch^
echo "Checking out the private branch..."
echo "[git checkout $PRIVATE_BRANCH]"
git checkout $PRIVATE_BRANCH

# Merge the public branch into the private branch to incorporate any changes
#echo "Merging the public branch into the private branch..."
#echo "[git merge github/main --no-ff -m 'Merge public changes']"
#git merge github/main --no-ff -m "Merge public changes" --allow-unrelated-histories

# Create a new branch for public commits
echo "Creating a new branch for public commits..."
echo "[git checkout -b $PUBLIC_BRANCH]"
git checkout -b $PUBLIC_BRANCH

# Squash all commits into a single commit
echo "Squashing all commits into a single commit..."
echo "[git reset $(git commit-tree HEAD^{tree} -m 'Release $(date +%Y-%m-%d)')]"
git reset $(git commit-tree HEAD^{tree} -m "Release $(date +%Y-%m-%d)")

# Force push the squashed commit to GitHub
echo "Force pushing the squashed commit to GitHub..."
echo "[git push github $PUBLIC_BRANCH:main --force]"
git push github $PUBLIC_BRANCH:main --force

# Checkout back to the private branch
echo "Checking out the private branch..."
echo "[git checkout $PRIVATE_BRANCH]"
git checkout $PRIVATE_BRANCH

# Delete the temporary public branch
echo "Deleting the temporary public branch..."
echo "[git branch -D $PUBLIC_BRANCH]"
git branch -D $PUBLIC_BRANCH