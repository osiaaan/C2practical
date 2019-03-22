patchfile=$1/src$2.patch
if [ -e $patchfile ] ;
then
  patch -fp0 < $1/src$2.patch
  make doc
else
  echo "Could not find $patchfile"
fi
