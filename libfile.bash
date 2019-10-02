#!/bin/bash
#
# This script has been developed at Scientific Computing Service of
# the National Biotechnology Center of the Spanish Research Council 
# (CNB-CSIC) and it is therefore the property of CNB-CSIC.
#
# (c) José R. Valverde, CNB-CSIC. 2018-2019
# Released under a LGPL or EU-LGPL license
#
#---------------------------------------------------------------------
#
#	This file contains function definitions for BASH
#
#	In order to use these functions you should source this file (e.g.
# . libfile.bash).
#
# File processing
# ---------------
#
# *recursively* get the last file pointed to by a symlink
function lastlink() { 
    [ ! -h "$1" ] && echo "$1" || (local link="$(expr "$(command ls -ld -- "$1")" : '.*-> \(.*\)$')"; cd $(dirname $1); lastlink "$link" | sed "s|^\([^/].*\)\$|$(dirname $1)/\1|"); 
}


# get the real absolute path name of a file
# uses lastlink above to identify the pointed-to file and then cd's to its
# directory to get the path and add it to the name
function realfile () {
   if [ -x `which readlink` ] ; then
        readlink -e $1
    elif [ -x `which realpath` ] ; then
        realpath $1
    else
        local target=`lastlink $1`
        local tname=${target##*/}
        local tdir=`dirname $target`
        echo $(cd $tdir ; pwd -P)"/$tname"
    fi
}



# a cheap substitute for realpath/readlink
# if realpath/readlink exits, it will use them, otherwise it will
# try its best to get the absolute path name without de-symlinking it!!!
function abspath() {
    x=`which readlink`
    if [ -x "$x" ] ; then readlink -f $1 ; return ; fi
    x=`which realpath`
    if [ -x "$x" ] ; then realpath $1 ; return ; fi
    # generate absolute path from relative path
    # $1     : relative filename
    # return : absolute path
    if [ -d "$1" ]; then
        # dir
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        # file
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    else
        #echo "/dev/null"
        echo " ---FILE-DOES-NOT-EXIST---"
    fi
}


# obtain relative path to go from source to target
# i.e. the path to use to make a link to target inside of source
function relpath() {
    # both $1 and $2 are absolute paths beginning with /
    # returns relative path to $2($target) from $1($source)
    if [ $# -eq 2 ] ; then
       source="$1"
       target="$2"
    elif [ $# -eq 1 ] ; then
       source=`pwd`
       target="$1"
    else
       echo "---relpath error, usage: relpath [source] target---"
       return
    fi
    
    # example use (ensure a symlink contains the relative path to 
    # a file from our current position instead of its absolute path):
    #    f=`abspath $file`
    #    ln -sf `relpath $f` ./target.lnk
    #    ln -sf `relpath . $f` .
     
#    if [ ! -e "$source" -o ! -e "$target" ] ; then
#        echo "---NO-SUCH-FILE---"
#        return
#    fi
    if [ ! -d "$source" ] ; then
        echo "---NO-SUCH-FILE-$source---"
        return
    else
        source=`abspath $source`
    fi
    if [ ! -e "$target" ] ; then
        echo "---NO-SUCH-FILE-$target---"
        return
    else
        target=`abspath $target`
    fi

    # example:
    #    f=`abspath $file`
    #    ln -s `relpath `pwd` $f` f.lnk

    common_part=$source # for now
    result="" # for now

    while [[ "${target#$common_part}" == "${target}" ]]; do
        # no match, means that candidate common part is not correct
        # go up one level (reduce common part)
        common_part="$(dirname $common_part)"
        # and record that we went back, with correct / handling
        if [[ -z $result ]]; then
            result=".."
        else
            result="../$result"
        fi
    done

    if [[ $common_part == "/" ]]; then
        # special case for root (no common path)
        result="$result/"
    fi

    # since we now have identified the common part,
    # compute the non-common part
    forward_part="${target#$common_part}"

    # and now stick all parts together
    if [[ -n $result ]] && [[ -n $forward_part ]]; then
        result="$result$forward_part"
    elif [[ -n $forward_part ]]; then
        # extra slash removal
        result="${forward_part:1}"
    fi

    echo $result
}

function posix_relpath () {
    [ $# -ge 1 ] && [ $# -le 2 ] || return 1
    current="${2:+"$1"}"
    target="${2:-"$1"}"
    [ "$target" != . ] || target=/
    target="/${target##/}"
    [ "$current" != . ] || current=/
    current="${current:="/"}"
    current="/${current##/}"
    appendix="${target##/}"
    relative=''
    while appendix="${target#"$current"/}"
        [ "$current" != '/' ] && [ "$appendix" = "$target" ]; do
        if [ "$current" = "$appendix" ]; then
            relative="${relative:-.}"
            echo "${relative#/}"
            return 0
        fi
        current="${current%/*}"
        relative="$relative${relative:+/}.."
    done
    relative="$relative${relative:+${appendix:+/}}${appendix#/}"
    echo "$relative"
}


# print file size in human readable form
function prsize() {
    local f sz TB GB MB KB
    for f in "$@" ; do
        #f=$1
        # derefernce links (-L)
        #sz=`stat -L -c%s "$f"`
        sz=`stat -c%s "$f"`
        # powers of 2
        TB=1099511627776
        GB=1073741824
        MB=1048576
        KB=1024
        # powers of 10
        #TB=1000000000000
        #GB=1000000000
        #MB=1000000
        #KB=1000
        # GB?
        if [ $sz -gt $TB ] ; then
	    #echo $((sz / TB)) GB $f
            sz=`echo "scale=1; $sz/$TB" | bc`
            echo ${sz}T  $f
        elif [ $sz -gt $GB ] ; then
	    #echo $((sz / GB)) GB $f
            sz=`echo "scale=1; $sz/$GB" | bc`
            echo ${sz}G $f
        elif [ $sz -gt $MB ] ; then
	    #echo $((sz / MB)) MB $f
            sz=`echo "scale=1; $sz/$MB" | bc`
            echo ${sz}M $f
        elif [ $sz -gt $KB ] ; then
	    #echo $((sz / KB)) KB $f
            sz=`echo "scale=1; $sz/$KB" | bc`
            echo ${sz}K $f
        else
            echo ${sz}B $f
        fi
    done
}

