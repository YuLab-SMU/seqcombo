#!/bin/bash
set -o errexit -o nounset
addToDrat(){
  PKG_REPO=$PWD

  cd ..; mkdir drat; cd drat

  ## Set up Repo parameters
  git init
  git config --global user.email "guangchuangyu@gmail.com"
  git config --global user.name "GuangchuangYu"
  git config --global push.default simple

  ## Get drat repo
  git remote add upstream "https://$GH_TOKEN@github.com/GITHUB_USERNAME/drat.git"
  git fetch upstream 2>err.txt


  Rscript -e "drat::insertPackage('$PKG_REPO/$PKG_TARBALL', \
    repodir = './docs', \
    commit='Travis update: build $TRAVIS_BUILD_NUMBER')"
  git push 2> /tmp/err.txt

}
addToDrat
