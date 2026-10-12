# echo "cleanup BesEvtGen BesEvtGen-00-04-08 in /workfs2/bes/cloudsail/workarea/boss708/Generator"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtBesEvtGentempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtBesEvtGentempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=BesEvtGen -version=BesEvtGen-00-04-08 -path=/workfs2/bes/cloudsail/workarea/boss708/Generator  $* >${cmtBesEvtGentempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=BesEvtGen -version=BesEvtGen-00-04-08 -path=/workfs2/bes/cloudsail/workarea/boss708/Generator  $* >${cmtBesEvtGentempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtBesEvtGentempfile}
  unset cmtBesEvtGentempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtBesEvtGentempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtBesEvtGentempfile}
unset cmtBesEvtGentempfile
exit $cmtcleanupstatus

