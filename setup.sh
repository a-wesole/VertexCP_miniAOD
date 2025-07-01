#this script is based on git-cms-merge-topic - there is an added line to only merge HeavyIonsAnalysis directory the rest are not needed
cd $CMSSW_BASE
mv src VCA

#cd ..
mkdir src
cd src

INITOPTIONS="--upstream-only"
git cms-init $INITOPTIONS
cd $CMSSW_BASE/src
BRANCH=forest_CMSSW_13_2_X
COMMAND_NAME=cms-merge-topic
TEMP_BRANCH_WORD=$(echo $COMMAND_NAME | cut -d'-' -f2)
TEMP_BRANCH=${TEMP_BRANCH_WORD}-attempt
git fetch . +HEAD:$TEMP_BRANCH
GITHUB_USER=CmsHI
REPOSITORY=https://github.com/CmsHI/cmssw.git
COMMIT=$BRANCH
FULL_BRANCH=$GITHUB_USER/$BRANCH
git fetch -n $REPOSITORY +$COMMIT:$FULL_BRANCH
CURRENT_BRANCH=`git rev-parse --abbrev-ref HEAD`
git checkout $TEMP_BRANCH
MERGE_BASE_BRANCH=$CURRENT_BRANCH
MERGE_BASE=`git merge-base $FULL_BRANCH $MERGE_BASE_BRANCH`
git cms-sparse-checkout $DEBUG_OPT $MERGE_BASE $FULL_BRANCH
git sparse-checkout set HeavyIonsAnalysis/EventAnalysis HeavyIonsAnalysis/TrackAnalysis
git read-tree -mu HEAD
git merge $NO_COMMIT $MERGE_STRATEGY $STRATEGY_OPTION --no-ff -m "Merged $BRANCH from repository $GITHUB_USER with cms-merge-topic" $GITHUB_USER/$BRANCH
git checkout $CURRENT_BRANCH
git read-tree -mu HEAD
git merge --ff $TEMP_BRANCH
git branch -D $TEMP_BRANCH
cd $CMSSW_BASE
mv VCA/* src
