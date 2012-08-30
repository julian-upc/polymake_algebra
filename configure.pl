#  Copyright (c) 1997-2012
#  Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Darmstadt, Germany)
#  http://www.polymake.org
#
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation; either version 2, or (at your option) any
#  later version: http://www.gnu.org/licenses/gpl.txt.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#-------------------------------------------------------------------------------
# $Project: polymake $$Id$

@make_vars=qw( CXXflags LDflags Libs );

sub allowed_options {
   my ($allowed_options, $allowed_with)=@_;
   @$allowed_with{ qw( singular ) }=();
}

sub usage {
   print STDERR "  --with-singular=PATH  installation path of libsingular\n";
}

sub proceed {
   my ($options)=@_;
   if (defined (my $singular_path=$$options{singular})) {
      my $singular_inc="$singular_path/include";
      my $singular_lib=Polymake::Configure::get_libdir($singular_path, "singular");
      if (-f "$singular_inc/libsingular.h" && -f "$singular_lib/libsingular.$Config::Config{dlext}") {
	 $CXXflags="-I$singular_inc";
	 $LDflags="-L$singular_lib -Wl,-rpath,$singular_lib";
      } else {
	 die "Invalid installation location of libsingular: header file libsingular.h and/or library libsingular.$Config::Config{dlext} not found\n";
      }

   } else {
      local $Polymake::Configure::Libs="-lsingular -ldl -lpthread $Polymake::Configure::Libs";
      my $error=Polymake::Configure::build_test_program(<<"---");
#include "libsingular.h"
#include <string>
#include <dlfcn.h>
int main() {
  Dl_info dli;
  if (!dladdr((void*)&siInit,&dli))
    return -1;
  char* cpath = omStrDup(dli.dli_fname);
  siInit(cpath);
  return 0;
}
---
      if ($?==0) {
	 $error=Polymake::Configure::run_test_program();
	 if ($?) {
	    die "Could not run a test program checking for libsingular.\n",
	        "The complete error log follows:\n$error\n\n",
		"Please investigate the reasons and fix the installation.\n";
	 }
      } else {
	 die "Could not compile a test program checking for libsingular.\n",
             "The most probable reasons are that the library is installed at a non-standard location,\n",
	     "is not configured to build a shared module, or missing at all.\n",
	     "The complete error log follows:\n$error\n\n",
	     "Please install the library and specify its location using --with-singular option, if needed.\n",
	     "Please remember to enable shared modules when configuring the libsingular!\n";
      }
   }

   $Libs="-lsingular";
}
