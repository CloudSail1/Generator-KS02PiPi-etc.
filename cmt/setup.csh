# echo "setup BesEvtGen BesEvtGen-00-04-08 in /workfs2/bes/cloudsail/workarea/boss708/Generator"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtBesEvtGentempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtBesEvtGentempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=BesEvtGen -version=BesEvtGen-00-04-08 -path=/workfs2/bes/cloudsail/workarea/boss708/Generator  -no_cleanup $* >${cmtBesEvtGentempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=BesEvtGen -version=BesEvtGen-00-04-08 -path=/workfs2/bes/cloudsail/workarea/boss708/Generator  -no_cleanup $* >${cmtBesEvtGentempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtBesEvtGentempfile}
  unset cmtBesEvtGentempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtBesEvtGentempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtBesEvtGentempfile}
unset cmtBesEvtGentempfile
exit $cmtsetupstatus

