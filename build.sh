if [ -d .svn ] ;
then
  echo "do not run the build.sh script in the svn directory -"
  echo "it is for the final build of the students version."
  echo "Use INSTALL instead"
  exit 1
fi

VERSION=2.3

if [ ! -d dune-common ] ;
then
  tar xvzf dune-common-$VERSION.tar.gz
  ln -s dune-common-$VERSION dune-common
fi
if [ ! -d dune-geometry ] ;
then
  tar xvzf dune-geometry-$VERSION.tar.gz
  ln -s dune-geometry-$VERSION dune-geometry
fi
if [ ! -d dune-grid ] ;
then
  tar xvzf dune-grid-$VERSION.tar.gz
  ln -s dune-grid-$VERSION dune-grid
fi
if [ ! -d dune-alugrid ] ;
then
  tar xvzf dune-alugrid-$VERSION.tar.gz
  ln -s dune-alugrid-$VERSION dune-alugrid
fi

if [ ! -d c2-practical ] ;
then
  tar xvzf c2-practical-1.tar.gz
  ln -s c2-practical-1 C2-practical
fi

sed "s@MAINPATH@$PWD@g" config.opts.in > config.opts
./dune-common/bin/dunecontrol --opts=config.opts all

cd C2-practical
make doc

