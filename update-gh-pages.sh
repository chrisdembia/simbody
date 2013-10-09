#!/bin/sh
# http://sleepycoders.blogspot.se/2013/03/sharing-travis-ci-generated-files.html
#if [ "$TRAVIS_PULL_REQUEST" == "false" ]; then
echo -e "Starting to update gh-pages.\n"
make doxygen
cp -R html $HOME/html
cd $HOME
git config --global user.email "travis@travis-ci.org"
git config --global user.name "Travis"
git clone --quiet --branch=gh-pages https://${GH_TOKEN}@github.com/chrisdembia/simbody.git gh-pages > /dev/null
cd gh-pages
git rm -rf .
cp -Rf $HOME/html/* .
git add -f .
git commit -m "Travis build $TRAVIS_BUILD_NUMBER pushed to gh-pages."
git push -fq origin gh-pages > /dev/null
echo -e "Done updating gh-pages.\n"
#fi

