#!/bin/bash
#coding :utf-8
set -e

if [[ "$(whoami)" == "root" ]]; then
echo Please do not run this as root, too much liability.
echo Only give your password for installation if prompted. Thanks!
exit
fi

tf="" # space-separated list of extracted folders/tmp files, not including java...

cleanup() {
tf="$tf ${sources}/*.delete"
tf="$tf $tmp"
if [[ "$tf" != "" ]]; then
set -- $tf
#      printf "cleaning up...\n"
for f in $@; do
if [ -e $f ]; then
#            printf "\trm -R $f\n"
#            sudo rm -R $f
if ! rm -Rf $f &> $log_file; then continue;fi
fi
done
fi
}

trap "cleanup" INT

### GLOBAL ###

tmp=/tmp/.cloudbio_ver_check
sources=/usr/local/hadoop/cloudbio_software
log_file=${sources}/cloudbio_install.log

make_checked=0
gcc_checked=0

fetch_curl="curl --location --max-redirs 3" #sourceforge involves redirect to mirror
fetch_wget=wget

# define defaults
fetch_agent=$fetch_wget
md5_agent="" # MacOSX has 'md5'

# likely stable
java_source=http://javadl.sun.com/webapps/download/AutoDL?BundleId=49015
java_size="21628961  jre-6u26-linux-i586.bin"
java_md5sum="9a8718922965deedd662439fdc3bc467  jre.bin"


# involves redirect from sourceforge, mac version of curl has trouble/can't seem to handle
bwa_source=http://sourceforge.net/projects/bio-bwa/files/bwa-0.5.9.tar.bz2
bwa_size="118705  bwa-0.5.9.tar.bz2"
bwa_md5sum="b60719d7f6bc49fd39cbc4e0b6bd787e  bwa-0.5.9.tar.bz2"
samtools_source=http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2
samtools_size="379306  samtools-0.1.18.tar.bz2"
samtools_md5sum="71dab132e21c0766f0de84c2371a9157  samtools-0.1.18.tar.bz2"


check_bash() {
# bash as shell
if [[ "$SHELL" != "/bin/bash" ]]; then
printf "shell ${SHELL##*/} not supported, exiting...\n" >&2
exit 1
fi
}

check_software() {
# checks if given software is installed in PATH
if ! which $1 &> /dev/null;then
printf "missing $1...\n"
return 1
fi
return 0
}

check_basic_software() {
# check for software needed to run
for c in sudo tar bzip2 gzip; do
if ! check_software $c;then
return 1
fi
return 0
done
}

define_md5_agent() {
if check_software md5sum; then
md5_agent=md5sum
else
if check_software md5; then
md5_agent=md5 # ur in a Mac
fi
fi
}


define_fetch_agent() {
if ! check_software wget; then
if ! check_software curl; then
printf "missing download utility such as wget or curl...\n" >&2
exit 1
else
fetch_agent=$fetch_curl
fi
else
fetch_agent=$fetch_wget
fi
}

verify_integrity() {
name=$1
fetched_filename=$2

if [ -e $fetched_filename ]; then
if [[ "$md5_agent" != "" ]];then
eval known_sum='${'$name'_md5sum}'
this_sum=$($md5_agent $fetched_filename)
#            echo comparing: $this_sum with $known_sum
if [[ "$this_sum" != "$known_sum" ]]; then
printf "\tfailed checksum\n"
return 1
#            else
#                printf "\tpassed checksum\n"
fi
else
#            return 1
printf "\twarning: cannot create checksum, comparing by size only.\n"
this_size=$(du -b $fetched_filename)
eval known_size='${'$name'_size}'
if [[ "${known_size%%[[:space:]]*}" != "${this_size%%[[:space:]]*}" ]] || [[ "${known_size##*[[:space:]]}" != "${this_size##*[[:space:]]}" ]]; then
printf "\tfailed size check\n"
return 1
#            else
#                printf "\tpassed size check\n"
fi
fi
return 0
fi
return 1
}

get_sudo_auth() {
if ! sudo -v; then
printf "Problem with sudo permission, cannot install.\n"
exit 2
fi
}


extract() {
name=$1
comp_type=$2

printf "\textracting $name...\n"
if [[ "$comp_type" == "-g" ]]; then
tar xzf $name # assumes same directory
else
if [[ "$comp_type" == "-b" ]]; then
tar xjf $name
else
printf "compression type $comp_type not found\n"
return 1
fi
fi
if [[ "$?" != "0" ]]; then
return 2
fi
return 0
}

get_and_decompress() {
name=$1
source_url=$2
comp_type=$3
source_file=${source_url##*/}
if [ ! -d $sources ]; then
mkdir $sources          # download to, extract from, etc
fi
cd $sources

if ! verify_integrity $name $source_file; then
if [ -e $source_file ]; then
# move to not overwrite and not file.N for existing with wget/curl
mv $source_file ${source_file}.delete
fi
printf "\tfetching $name from $source_url...\n"
if [[ "$fetch_agent" == "$fetch_curl" ]];then
fetch_command="$fetch_curl --remote-name"
else
fetch_command=$fetch_agent
fi
if ! $fetch_command $source_url &> $log_file; then
return 1
fi
if ! verify_integrity $name $source_file; then
return 1
fi
else
printf "\tfound intact $source_file in $(pwd)\n"
prev_extracted_dir=$(ls ${source_file:0:3}* | head -1)
if [ -d $prev_extracted_dir ]; then
rm -R $prev_extracted_dir
fi
fi

if ! extract $source_file $comp_type; then
return 1
fi
source_dir=$(ls -d ${source_file:0:3}* | head -1) # dir usually shorter than tarball
tf="$tf ${sources}/${source_dir}"
if [ ! -d $source_dir ]; then
return 1
fi
if [[ ${source_dir:0:3} != ${source_file:0:3} ]]; then
printf "cannot find match for extracted directory.\n" >&2
cd - >& /dev/null
return 1
fi
cd $source_dir
return 0
}

### software_check.sh ###


check_bwa() {
printf "BWA\t\t"
if which bwa &> /dev/null; then
bwa 2> $tmp
bwa_ver=$(awk '{if ($1 == "Version:") {print substr($2,3,3);exit} }' $tmp)
if [[ "${bwa_ver%.*}" -ge "5" && "${bwa_ver#*.}" -ge "7" ]]; then
printf "installed.\n"
return 0
else

printf "bad version.\n"
return 1
fi
else
printf "not installed.\n"
return 2
fi
}

check_samtools() {
printf "samtools\t"
if which samtools &> /dev/null; then
samtools &> $tmp
samtools_ver=$(awk '{if ($1 == "Version:") {print substr($2,3,4);exit} }' $tmp)
if [[ "${samtools_ver%.*}" -ge "1" && "${samtools_ver#*.}" -ge "10" ]]; then
printf "installed.\n"
return 0
else
printf "bad version.\n"
return 1
fi
else
printf "not installed.\n"
return 2
fi
}

check_java() {
printf "Java\t\t"
if java -version &> /dev/null; then
java -version >& $tmp
java_ver=$(awk '{if ($1 == "java") {print substr($3,2,3);exit} }' $tmp)
if [[ "${java_ver%.*}" -ge "1" && "${java_ver#*.}" -ge "5" ]]; then
printf "installed\n"
return 0
else
printf "bad version\n"
return 1
fi
else
printf "not installed.\n"
return 2
fi
}


#    get_sudo_auth


do_java() {
if ! check_java; then
mkdir -p "$sources"
cd "$sources"
if ! verify_integrity java jre.bin; then
if [[ "$fetch_agent" == "$fetch_curl" ]]; then
fetch_command="$fetch_curl -o jre.bin $java_source";
else
fetch_command="$fetch_wget -O jre.bin $java_source";
fi
printf "\tfetching java from $java_source...\n"
if ! $fetch_command &> $log_file; then
printf "Problem fetching jre.bin. Aborting...\n"
return 1
fi
fi
cd /usr/local/bin
printf "\tbuilding java...\n"
# as per java's instructions, go to dir to install, call script from there.
if ! sudo sh ${sources}/jre.bin 1>> $log_file; then
return 1
fi
java_dir=$(ls -d /usr/local/bin/jre* | head -1)
if [ ! -d $java_dir ]; then
return 1
fi
# this still would require an update to PATH, so just creating symlink from a spot already in PATH
printf "\tcreating symbolic link from ${java_dir}/bin/java to /usr/local/bin/java\n"
sudo ln -s "${java_dir}/bin/java" "/usr/local/bin/java" # should be '.', but want to be explicit
cd - >& /dev/null
fi
# printf "done.\n"
return 0
}

do_bwa() {
if ! check_bwa; then
if [ $make_checked -ne 1 ]; then
if ! check_software make; then
return 1
else
make_checked=1
fi
fi
if [ $gcc_checked -ne 1 ]; then
if ! check_software gcc; then
return 1
else
gcc_checked=1
fi
fi
if ! get_and_decompress "bwa" $bwa_source -b; then
return 1
fi

#    get_sudo_auth

# now dropped into newly extracted directory!
printf "\tbuilding bwa...\n"
if ! make &> $log_file; then
return 1
fi
printf "\tcopying bwa binary to /usr/local/bin...\n"
if ! sudo cp bwa /usr/local/bin; then
return 1
fi
fi
# printf "done.\n"
return 0
}

do_samtools() {
if ! check_samtools; then
if [ $make_checked -ne 1 ]; then
if ! check_software make; then
return 1
else
make_checked=1
fi
fi
if [ $gcc_checked -ne 1 ]; then
if ! check_software gcc; then
return 1
else
gcc_checked=1
fi
fi

if ! get_and_decompress "samtools" $samtools_source -b; then
return 1
fi

# now dropped into newly extracted directory!
printf "\tbuilding samtools...\n"
if ! make &> $log_file; then
return 1
fi
printf "\tcopying samtools binary to /usr/local/bin...\n"
if ! sudo cp samtools misc/*.pl /usr/local/bin; then
return 1
fi
fi
#printf "done.\n"
return 0
}


# list only those needed to be installed..., add to growing list
do_checks() {
echo "Checking for software needed for Cloudbio..."
printf "\t\tStatus\n"
for software in bwa samtools java; do
if ! check_${software};then continue;fi;
done
}

# order is important!
# python needed by python modules, samtools for pysam for other python modules, etc
do_install() {
echo "Installing software needed for Cloudbio..."
printf "\t\tStatus\n"
for software in bwa samtools python py_packs java; do
if ! do_${software}; then
echo Something untoward occured with $software
echo Please check $log_file for details.
cleanup
exit -1
fi
done
cleanup
}

check_bash           # may be somewhat redundant...
define_md5_agent     # MacOSX has 'md5', not 'md5sum'
define_fetch_agent   # MacOSX has curl, not wget
check_basic_software # do we have tar, bzip2, gzip and sudo?
# define_python_name occurs with check_python for python or python2

case $1 in
-i)
do_install
;;
-c)
do_checks
;;
*)
echo "software check
Usage: $0 [OPTION]
Determine system state to run DNA_Seq, installing software if needed.
-h  print this help message
-i  check for and install software
(prints to $log_file)
-c  check software only"
;;
esac
