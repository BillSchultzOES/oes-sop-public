#!/bin/sh

set -e

#git config --global user.email "jake@jakebowers.org"
#git config --global user.name "Jake Bowers"

cd ..
rm -rf book-output
git clone -b gh-pages git@github.com:gsa-oes/sop.git book-output
cd book-output
cp -r ../Book/_book/* ./
git add --all *
git commit -m "Update the book" || true
git push -q origin gh-pages
git checkout main 
