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
echo "Merging the public branch into the private branch..."
echo "[git merge github/main --no-ff -m 'Merge public changes']"
git merge github/main --no-ff -m "Merge public changes" --allow-unrelated-histories
