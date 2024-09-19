#! /bin/sh
func() {
    echo "Usage:"
    echo "single_fa.sh [-m mode(many or only)] [-n header name (if -m=only)] [-i input_file] [-o output_file]"
    echo "MODE:"
    echo "many: keep all header"
    echo "only: keep single header"
    exit -1
}

while getopts :m:n:i:o: varname
do
   let optnum++
   case $varname in
   m)
      mode="$OPTARG" ;;
   n)
      name="$OPTARG" ;;
   i)
      input="$OPTARG" ;;
   o)
      output="$OPTARG" ;;
   ?) 
      func ;;
   esac
done

if [ $mode = "many" ];then
     awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'  $input > $output
elif test  "$name" ; then
     sed '1!{/^>/d;}' $input | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' |sed  "1i\>$name" > $output
else 
     sed '1!{/^>/d;}' $input | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'  > $output
fi
