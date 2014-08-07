#!/usr/bin/perl
#use warnings;        # Avertissement des messages d'erreurs
#use strict;          # Verification des declarations
use File::Spec::Functions;
use File::Basename;
use URI::file;
use Cwd "realpath";
use Getopt::Long;

#########################################################################
#                                                                       #
#                      List of global variables                         #
#                                                                       #
#########################################################################

my $msg =  2;
my %msgl = ("q", 0 , "e", 1 ,"w", 2 ,"v", 3 ,"vv", 4,"vvv", 5 ) ;
my $items_per_line = 4 ;   # number of items per Makefile line
my $item_count = 0;
my $ext = undef;
my @listfile;
my @includelist;
my $use_strict = undef;
my $deep_include = undef;
my $soft_restriction = undef;
my $flat_layout = undef;
my $short_target_names = undef;
my $local_dir = undef;
my $dup_ok = undef;
my $side_dir_inc = undef;
my $anywhere_inc = undef;
my $export_list = undef;
my @current_dependencies_list = ();
my @outside_tree_list = ();
my %module_missing = ();
my %module_missing_ignored = ();
my %include_missing_ignored = ();
my $bool = 0;
my %LISTOBJECT = ( ); # Hash of SRCFile object with filename, path and extension as the key

#########################################################################
#                                                                       #
#              List of function and object definition                   #
#                                                                       #
#########################################################################

{ package SRCFile;
    {
    
    # List of the type of file associated with the extension
    our %TYPE_LIST = (
        f     => "COMPILABLE",
        ftn   => "COMPILABLE",
        ptn   => "INCLUDE",
        f90   => "COMPILABLE",
        ftn90 => "COMPILABLE",
        ptn90 => "INCLUDE",
        cdk   => "INCLUDE",
        cdk90 => "COMPILABLE",
        c     => "COMPILABLE",
        h     => "INCLUDE",
        hf    => "INCLUDE",
        fh    => "INCLUDE",
        inc   => "INCLUDE",
        tmpl90 => "COMPILABLE",
    );

    # @new: Constructor of the SRCFile Object
    # 
    # input:
    #   $1 = {
    #           path => 
    #           filename =>
    #           extension =>
    #        }
    #
    # output:
    #   pointer to the object

    sub new 
    {
        my ( $class, $ref_arguments ) = @_;
        
        $class = ref($class) || $class;
        my $this = {};
        bless( $this, $class );
            
        $this->{FULLPATH_SRC}     = " ";
        $this->{FILENAME}         = " ";
        $this->{EXTENSION}        = " ";

        $this->{FULLPATH_SRC}     = $ref_arguments->{path}; 
        $this->{FILENAME}         = $ref_arguments->{filename};
        $this->{EXTENSION}        = $ref_arguments->{extension};
      %{$this->{DEPENDENCIES}}    = ();   # A list to the required file
        $this->{TYPE}             = $TYPE_LIST{lc $this->{EXTENSION}};   # Type based on the extension
        $this->{STATUS}           = undef;
      @{$this->{UNSOLVED_MODULE}} = ();
      @{$this->{MODULE_LIST}}     = ();
      %{$this->{UNKNOWN_MODULE}}  = ();
      %{$this->{UNKOWN_USE}}      = ();
            
        return $this;
    }

    # @displ: Display file information
    #
    # input:
    #   none
    #
    # output:
    #   none
    sub displ {
        my $self = $_[0];
        my $key;
        print "\t$_[0]->{FULLPATH_SRC} \t$_[0]->{FILENAME} \t$_[0]->{EXTENSION} \t$_[0]->{TYPE}\t\t$_[0]->{STATUS}\n";
        for (keys %{$self->{DEPENDENCIES}}) {
            print "\t\tDEPENDENCIES:\t $_\n";
        }
        for (@{$self->{MODULE_LIST}}) {
            print "\t\tMODULE:\t $_\n";
        }
        for (@{$self->{UNSOLVED_MODULE}}) {
            print "\t\tMISSING MODULE:\t $_\n";
        }
    } 

    # @getFilename: Get full filename with path and extension
    # 
    # input: 
    #   none
    #
    # output:
    #   full path
    sub getFilename {
        return "$_[0]->{FULLPATH_SRC}$_[0]->{FILENAME}.$_[0]->{EXTENSION}";
    }

    # @set_status: Set status of object to true
    #
    # input: none
    #
    # output: none
    sub set_status {
        $_[0]->{STATUS} = 1; #true
    }

    # @get_status: get status the object
    #
    # input: none
    #
    # output: 
    #   status
    sub get_status {
        return $_[0]->{STATUS};
    }

    # @reset_status: Set status of object to false
    #
    # input: none
    #
    # output: none
    sub reset_status {
        $_[0]->{STATUS} = undef; #false
    }

    # @has_module: find if the object defined the module 
    #
    # intput: 
    #   $1 = Module name to find
    #
    # output:
    #   true (1) if the module as been found, false (undef) otherwise
    sub has_module {
        my $module_name = $_[1];    
        for (@{$_[0]->{MODULE_LIST}}) {
            return 1 if $_ eq $module_name;
        }
        return undef;
    }

    # @has_unsolved_module: find if the object has the module in his unsolved module list
    #
    # input: 
    #   $1 = Module name to find
    #
    # output:
    #   true (1) if the module as been found, false (undef) otherwise
    sub has_unsolved_module {
        my $module_name = $_[1];        
        for (@{$_[0]->{UNSOLVED_MODULE}}) {
            return 1 if ($_ eq $module_name);
        }
        return undef;
    }

    # @remove_unsolved_module: delete the module in the unsolved module list
    #
    # input:
    #   $1 = Module name to delete
    #
    # output: none
    sub remove_unsolved_module {
        my $module_name = "$_[1]";
        my $idx = 0;
        for my $module (@{$_[0]->{UNSOLVED_MODULE}}) {
            if ($module eq $module_name) { 
                delete ${$_[0]->{UNSOLVED_MODULE}}[$idx]; 
            } 
            $idx++; 
        }
    }

    # @find_depedencies: find if the object has the filename in his depedencie list
    #
    # input:
    #   $1 = filename to search
    # 
    # output
    #   true (1) if the module as been found, false (undef) otherwise
    sub find_depedencies {
        my $search_depedencies = "$_[1]";
        #while(my($dep_filename, $dep_ptr) = each(%{$_[0]->{DEPENDENCIES}})) {
		  for my $dep_filename (keys %{$_[0]->{DEPENDENCIES}}) {
				my $dep_ptr = ${$_[0]->{DEPENDENCIES}}{$dep_filename};
            return 1 if (($search_depedencies eq $dep_filename) and ($_[0]->{FULLPATH_SRC} ne $dep_ptr->{FULLPATH_SRC} ) );
            return undef if (($search_depedencies eq $dep_filename) and ($_[0]->{FULLPATH_SRC} eq $dep_ptr->{FULLPATH_SRC} ) );
            return 1 if ($dep_ptr->SRCFile::find_depedencies($search_depedencies) );
        }
        return undef;
    }

    } #end package SRCFile
}


# @reset_all_file: Reset status of all object
#
# input:
#   %0 = Hash of object
#
# output: none
sub reset_all_file {
	 for (keys %{$_[0]}) {
		  $LISTOBJECT{$_}->reset_status();
	 }
}

# @search_undone_file: Find the key of the object that hasn't been processed yet (based on its status)
#
# input: 
#   %0 = Hash of object
#
# output:
#   key of the object, undef otherwise.
sub search_undone_file {
	 #$_[0] is ignored
	 for (keys %LISTOBJECT) {
		  return $_ if !$LISTOBJECT{$_}->get_status(); 
	 }
	 return undef;
}

# @search_module: Find the key of the object that own the module name
#
# input: 
#   %0 = Hash of objects
#   $1 = Module name
#
# output:
#   key of the object where the module has been found, undef otherwise.
sub search_module {
	 #$_[0] is ignored
    my $module_name = $_[1];
	 for (keys %LISTOBJECT) {
		  return $_ if ($LISTOBJECT{$_}->has_module($module_name)); 
    }
    return undef;
}

# @search_unsolved_module: Find the key of the first object that has the module as one of his unsolved module list
#
# input:
#   %0 = Hash of object
#   $1 = Module name to be search
#
# output:
#   key of the object, undef otherwise.
sub search_unsolved_module {
    my $module_name = $_[1];
	 for (keys %LISTOBJECT) {
		  return $_ if ($LISTOBJECT{$_}->SRCFile::has_unsolved_module($module_name));
	 }
	 return undef;
}

# @print_header: Print the first line of a dependency rule or variable list
#
# input:
#   $0 = First word(s) of line
#   $1 = Seperator
#   $2 = Word/item right after seperator, empty if no word is needed
#
# output: none
sub print_header {
	 $item_count = 0 ;
	 my($item1,$separ,$item2) = @_ ;
	 print STDOUT "$item1$separ" ;
	 print STDOUT "\t$item2" if ( "$item1" ne "$item2" && "$item2" ne "" ) ;
}

# @print_item: print each item of dependency rule or variable list (items_per_line items per line)
#
# input:
#   $0 = Item to print
#
# output: none
sub print_item {
	 my($item1) = @_ ;
	 if ($item1) {
		  print STDOUT " \\\n\t" if ($item_count == 0) ;
		  print STDOUT "$item1  ";
		  $item_count = 0 if ($item_count++ >= $items_per_line);
	 }
}

# @find_same_filename: Look if the filename is already used somewhere else.
#
# input: 
#   %0 = Hash of objects
#   $1 = filename to compare with.
#
# output:
#    key (filename) of the object if the file already exist in the list, false (undef) otherwise.
sub find_same_filename {
	 #$_[0] is ignored
    my $cmp_file = $LISTOBJECT{$_[1]};
    # Return if the soft_restriction is activated and if the file isn't compilable. 
    # soft_restriction = disable warning if 2 headers files have the same name
    return undef if ($soft_restriction and ($cmp_file->{TYPE} ne "COMPILABLE"));

	 for (keys %LISTOBJECT) {
        return $_ if (
				("$LISTOBJECT{$_}->{FILENAME}.$LISTOBJECT{$_}->{EXTENSION}" 
				 eq "$cmp_file->{FILENAME}.$cmp_file->{EXTENSION}") 
				and ($LISTOBJECT{$_}->{FULLPATH_SRC} ne $cmp_file->{FULLPATH_SRC}));
	 }
    return undef;
}

# @find_same_filename: Look if the filename is already used somewhere else.
#
# input: 
#   %0 = Hash of objects
#   $1 = path to filename to compare with.
#   $2 = filename.ext to compare with.
#
# output:
#    key (filename) of the object if the file already exist in the list, false (undef) otherwise.
sub find_same_filename2 {
	 #$_[0] is ignored
    my $mydir = $_[1];
	 my $myfilename = $_[2];
	 for my $key (keys %LISTOBJECT) {
		  my $file = $LISTOBJECT{$key};
		  if (("$file->{FILENAME}.$file->{EXTENSION}" eq $myfilename) 
				and ($file->{FULLPATH_SRC} ne $mydir)) {
				# Return if the soft_restriction is activated and if the file isn't compilable. 
				# soft_restriction = disable warning if 2 headers files have the same name
				return undef if ($soft_restriction and ($file->{TYPE} ne "COMPILABLE"));
				return $key
		  }
    }
    return undef;
}

# @find_same_output: Look if the filename is already used in the Object list.
#
# input: 
#   %0 = Hash of objects
#   $1 = Object to compare with.
#
# output:
#   key (filename) of the object if the file already exist in the list, false (undef) otherwise.
sub find_same_output {
	 #$_[0] is ignored
    my $cmp_file = $LISTOBJECT{$_[1]};
    return undef if $cmp_file->{TYPE} ne "COMPILABLE";
	 for my $key (keys %LISTOBJECT) {
        return $key if (
				($LISTOBJECT{$key}->{FILENAME} eq $cmp_file->{FILENAME}) 
				and ($LISTOBJECT{$key}->{TYPE} eq "COMPILABLE") 
				and ($key ne $_[1])); 
    }
    return undef;
}

# @find_string_in_array:
# 
# input: 
#   @0 = Array to search in.
#   $1 = String to search.
#
# output:
#   1 if the string as been found, undef otherwise.
sub find_string_in_array {
    my @myArray = @{$_[0]}; 
    my $string = $_[1];
    for my $string_tmp (@myArray) {
        return 1 if ($string eq $string_tmp);
    }
    return undef;
}

# @rec_print_dependencies: 
#
# input:
#   %0 = Hash of objects
#   $1 = Filename to print dependencies
#
# output:
#   none
sub rec_print_dependencies {
    my $file = ${$_[0]}{$_[1]};
    
	 #print STDERR "rec_print_dependencies: $file->{FILENAME} : $_[1]\n" if ($msg >= 5);
    #while(my($dep_filename, $dep_ptr) = each(%{$file->{DEPENDENCIES}})) {
	 for my $dep_filename (sort keys %{$file->{DEPENDENCIES}}) {
	 	  my $dep_ptr = ${$file->{DEPENDENCIES}}{$dep_filename};

		  my $tmp_filename = $dep_filename;
        $tmp_filename = "$dep_ptr->{FULLPATH_SRC}$dep_ptr->{FILENAME}.o" if ($dep_ptr->{TYPE} eq "COMPILABLE");
		  my $tmp_filename0 = $tmp_filename;
		  if ($flat_layout) {
				$tmp_filename0 = "$dep_ptr->{FILENAME}.$dep_ptr->{EXTENSION}";
				$tmp_filename0 = "$dep_ptr->{FILENAME}.o" if ($dep_ptr->{TYPE} eq "COMPILABLE");
		  }
    
		  #print STDERR "$file->{FILENAME}: $dep_filename : $_[1] : $tmp_filename0 : $tmp_filename\n" if ($msg >= 5);
        next if (($_[1] eq $dep_filename) or find_string_in_array(\@current_dependencies_list, $tmp_filename) );
        
        print_item($tmp_filename0);
        push @current_dependencies_list, $tmp_filename;

        # Recursively call the function to print all depedencies
        rec_print_dependencies(\%{$_[0]}, $dep_filename) if ($dep_ptr->{TYPE} ne "COMPILABLE" );
    }
}

# @has_legal_extension: 
#
# input: 
#   $0 = Extension to search
#
# output:
#   1 if the extension is valid, undef otherwise.
sub has_legal_extension {
    my $search_extension = lc  $_[0];
    for my $extension (keys(%SRCFile::TYPE_LIST)) {
        return 1 if $extension eq $search_extension;
    }
    return undef;
}


# @process_file
#
# input: 
#   $0 = file
#   $1 = ==1 if duplicatedfile_ok
#
# output:
#   undef if ok; 1 otherwise
sub pre_process_file {
	 my $entry = $_[0];
	 my $_dup_ok = $_[1];
	 return 1 if (! -f "$entry");

	 my $file = "$entry" ;
	 $file =~ s/,v$// ;
	 $file =~ s/[\s]+// ;
	 $file = File::Spec->abs2rel( canonpath($file), "./") if $file =~ /^\//;
    
	 return 1 if $file !~  /(.*\/)*(.*)[.]([^.]*$)/ ;  # ,v and path trimmed filename must be optional/path/root_file.extension
	 return 1 if exists $LISTOBJECT{$file}; 
	 
	 my $path = ($1 ? $1 : "");
	 my $filn = ($2 ? $2 : "");
	 my $exte = ($3 ? $3 : "");
	 # my $duplicated_filename = "";
	 
	 return 1 if !has_legal_extension($exte);

	 my $duplicated_filename1 = find_same_filename2(\%LISTOBJECT, $path, "$filn.$exte");

	 if ($duplicated_filename1 and $_dup_ok) {
		  delete $LISTOBJECT{$duplicated_filename1};
	 }

	 $LISTOBJECT{"$path$filn.$exte"} = new SRCFile({path => $path, filename => $filn, extension => $exte});

	 if ($duplicated_filename1 and $_dup_ok) {
		  print STDERR "WARNING: $duplicated_filename1 was replaced by $path$filn.$exte : ".$LISTOBJECT{"$path$filn.$exte"}->{FILENAME}.$LISTOBJECT{"$path$filn.$exte"}->{STATUS};
	 }
    
	 # Error handler
	 my $duplicated_filename2 = find_same_output(\%LISTOBJECT, "$path$filn.$exte");
	 if ($_dup_ok) {
		  if ($msg >= 1) {
				print STDERR "WARNING: using 2 files with the same name $duplicated_filename1 with $path$filn.$exte\n" if ($duplicated_filename1);
				print STDERR  "WARNING: using 2 files ($duplicated_filename2 and $path$filn.$exte) that will produce the same object file ($filn.o)\n" if ($duplicated_filename2);
		  }
	 } else {
		  die "ERR: using 2 files with the same name $duplicated_filename1 with $path$filn.$exte" if ($duplicated_filename1);
		  die "ERR: using 2 files ($duplicated_filename2 and $path$filn.$exte) that will produce the same object file ($filn.o)\n" if ($duplicated_filename2);		  
	 }

	 # print STDERR "process: '$entry' dupok=$_dup_ok ; path=$path ; filen=$filn ; exte=$exte ; dup=$duplicated_filename1\n" if ($msg >= 5);
	 
	 return undef;
}


# @find_inc_file
#
# input: 
#   $0 = file obj
#   $1 = supposed path to file to include
#   $2 = filename to include
#
# output:
#   actual path to file; undef if not found
sub  find_inc_file {
	 my $myfile = $_[0];
	 my $mypath = $_[1];
	 my $myfilename = $_[2];
	 for (@includelist) {
		  if (-f $_.'/'.$mypath.$myfilename) {
				return $_.'/'.$mypath;
		  } elsif (-f $_.'/'.$myfilename) {
				return $_.'/';
		  }
	 }
	 if ($anywhere_inc) {
		  for (keys %LISTOBJECT) {
				my $myobj=$LISTOBJECT{$_};
				my $myfilename2 = "$myobj->{FILENAME}.$myobj->{EXTENSION}";
				if ($myfilename eq $myfilename2 and -f "$myobj->{FULLPATH_SRC}$myfilename2") {
					 return $myobj->{FULLPATH_SRC};
				}
		  }
	 }
	 if ($side_dir_inc) {
		  my @mydirs = File::Spec->splitdir($myfile->{FULLPATH_SRC});
		  for my $mysubdir ('*','*/*','*/*/*','*/*/*/*','*/*/*/*/*') {
				my @myfile2 = glob "$mydirs[0]/$mysubdir/$myfilename\n";
				if ($myfile2[0]) {
					 # print STDERR "Found $myfilename in ".dirname($myfile2[0]) if ($msg >=4);
					 return dirname($myfile2[0]).'/';
				}
		  }
	 }
	 return undef;
}


# @process_file_for_include
#
# input: 
#   $0 = file object
#   $1 = filename
#
# output:
#   undef if ok; 1 otherwise
sub process_file_for_include {
	 my $file = $_[0];
	 my $tmp_dir = $_[1];
	 my $include_path = "";

	 if ($tmp_dir =~ /^\.\.\//) {
		  $include_path = File::Spec->abs2rel( canonpath("$file->{FULLPATH_SRC}/$tmp_dir"), "./"); # Convert file path position relatively to the base path
	 } elsif (-f canonpath("$file->{FULLPATH_SRC}/$tmp_dir")) {
		  $include_path = File::Spec->abs2rel( canonpath("$file->{FULLPATH_SRC}/$tmp_dir"), "./");
	 } else {
		  $include_path = File::Spec->abs2rel( canonpath($tmp_dir), "./");
	 }

	 # print STDERR "Missing $file->{FILENAME}.$file->{EXTENSION}: $tmp_dir\n" if (!$include_path and $msg>=4);

    # ,v and path trimmed filename must be optional/path/root_file.extension
	 
	 if ($include_path !~  /(.*\/)*(.*)[.]([^.]*$)/) {
		  # print STDERR "Outside $file->{FILENAME}.$file->{EXTENSION}: $tmp_dir : $include_path\n" if ($msg>=4);
		  return 1;
	 }

	 my $path = ($1 ? $1 : "");
	 my $filn = ($2 ? $2 : "");
	 my $exte = ($3 ? $3 : "");
	 my $duplicated_filename = "";

	 if (!has_legal_extension($exte)) {
		  # print STDERR "Bad Extension $file->{FILENAME}.$file->{EXTENSION}: $tmp_dir : $include_path : $exte\n" if ($msg>=4);
		  return 1;
	 }
	 if (! -f "$path$filn.$exte") {
		  if (!find_string_in_array(\@outside_tree_list, "$path$filn.$exte")) {
				my $path1 = find_inc_file($file,$path,"$filn.$exte");
				if (!$path1) {
					 # print STDERR "No file $file->{FILENAME}.$file->{EXTENSION}: $tmp_dir : $include_path : $path$filn.$exte\n" if ($msg>=4);
					 if (!exists($include_missing_ignored{"$path$filn.$exte"})) {
						  push @outside_tree_list, "$path$filn.$exte" if (!find_string_in_array(\@outside_tree_list, "$path$filn.$exte"));
					 }
					 return 1;
				}
				# print STDERR "Found $filn.$exte in $path1\n" if ($msg >=5);
				$path = $path1;
		  } else {
				return 1;
		  }
	 }

	 # Add file in the database if it's not in yet and if the file really exists.
	 $LISTOBJECT{"$path$filn.$exte"} = new SRCFile({path => $path, filename => $filn, extension => $exte}) 
		  if (!exists $LISTOBJECT{"$path$filn.$exte"});

	 # Force the file to not be analysed.
	 $LISTOBJECT{"$path$filn.$exte"}->set_status() 
	 	  if (!$deep_include);

	 # Error handler
	 die "ERR: using 2 files with the same name $duplicated_filename with $path$filn.$exte\n" 
		  if ($duplicated_filename = find_same_filename(\%LISTOBJECT, "$path$filn.$exte"));
	 die "ERR: using 2 files ($duplicated_filename and $path$filn.$exte) that will produce the same object file ($filn.o)\n" 
		  if ($duplicated_filename = find_same_output(\%LISTOBJECT, "$path$filn.$exte"));
	 die "ERR: cannot include compilable file ($tmp_dir) in $tmp_dir while using strict mode\n" 
		  if ($use_strict and $LISTOBJECT{"$path$filn.$exte"}->{TYPE} eq "COMPILABLE");

	 # Add to dependencies, if not already there
	 ${$file->{ DEPENDENCIES }}{"$path$filn.$exte"} = $LISTOBJECT{"$path$filn.$exte"} if (!exists ${$file->{ DEPENDENCIES }}{"$path$filn.$exte"});

	 return undef;
}


#########################################################################
#                                                                       #
#                        Main program beginning                         #
#                                                                       #
#########################################################################

#
# Process command line arguments
#
my $help = 0;
my $output_file='';
my $include_dirs='';
my $suppress_errors_file='';
GetOptions('verbose:+' => \$msg,
			  'quiet' => sub{$msg=0;},
			  'help' => \$help,
			  'flat_layout' => \$flat_layout,
			  'short_target_names' => \$short_target_names,
			  # 'local' => \$local_dir,
			  # 'dup_ok' => \$dup_ok,
			  'side_dir_inc' => \$side_dir_inc,
			  'any_inc' => \$anywhere_inc,
			  'strict' => \$use_strict,
			  'deep-include' => \$deep_include,
			  'soft-restriction' => \$soft_restriction,
			  'exp=s' => \$export_list,
			  'out=s' => \$output_file,
			  'includes=s' => \$include_dirs,
           'supp=s' => \$suppress_errors_file,
	 )
	 or $help=1;

@listfile = (@ARGV);
if (!$help and !($#listfile+1)) {
	 $help=1;
	 print STDERR "ERRROR: you must provide a list of targets\n"
}
if ($help) {
	 print STDERR "
Usage: s.dependencies.pl [-v|--quiet] \\
                         [--strict] [--deep-include] [--soft-restriction] \\
                         [--flat_layout] [--short_target_names] \\
                         [--exp=output_of_produced_file] [--out=outfile] \\
                         [--side_dir_inc] [--any_inc] \\
                         [--includes=list_of_inc_dirs]  \\ 
                         list_of_targets
       list_of_targets : must be a list of files or dirs
       list_of_inc_dirs: must be a list of ':'-separated dirs\n\n";
	 exit;
}

print STDERR "
s.dependencies.pl \\
   -v=$msg --strict=$use_strict --deep-include=$deep_include --soft-restriction=$soft_restriction \\
   --flat_layout=$flat_layout --short_target_names=$short_target_names \\
   --exp=$export_list --out=$output_file \\
   --side_dir_inc=$side_dir_inc --any_inc=$anywhere_inc \\
   --includes=$include_dirs \\
   ...
   \n" if ($msg>=3);

if ($output_file) {
	 print STDERR "Redirecting STDOUT to $output_file\n" if ($msg>=3);
	 open(STDOUT,">", "$output_file") or die "ERROR: Can't redirect STDOUT\n";
}
@includelist = split(':',$include_dirs) if ($include_dirs);

#
# Pre-Processsuppress_errors_file
#
if ($suppress_errors_file) {
	 print STDERR "process suppress_errors_file: $suppress_errors_file\n" if ($msg >= 3);
	 if (!open(INPUT,"<", $suppress_errors_file)) {
		  print STDERR "ERROR: Can't open supp file, ignoring: ".$suppress_errors_file."\n";
	 } else {
		  while (<INPUT>) {
				if ($_ =~ /^[\s]*module_missing[\s]+([^\s]+)/i) {
					 print STDERR "Suppressing missing mod msg for: ".$1."\n" ;#if ($msg >= 3);
					 $module_missing_ignored{$1} = 1;
				} elsif ($_ =~ /^[\s]*include_missing[\s]+([^\s]+)/i) {
					 print STDERR "Suppressing missing inc msg for: ".$1."\n" ;#if ($msg >= 3);
					 $include_missing_ignored{$1} = 1;
				} else {
					 print STDERR "Ignoring supp file line: ".$_."\n" ;#if ($msg >= 3);
				}
		  }
		  close INPUT;
	 }
}

#
# Pre-Process/register files
#
for (@listfile){
	 if (-d $_) {
		  print STDERR "process: '$_' $dup_ok\n" if ($msg >= 3);
		  for (glob "$_/*") {
				pre_process_file($_,$dup_ok);
		  }
	 }
	 else {
		  print STDERR "process: '$_'\n" if ($msg >= 3);
		  for (glob "$_") {
				pre_process_file($_,$dup_ok);
		  }
	 }
}

$dup_ok = 1;
if ($local_dir) {
	 for (glob "./*") {
		  pre_process_file($_,$dup_ok);
	 }
}

#
# Process files
#
reset_all_file(\%LISTOBJECT);

while(my $filename = search_undone_file(\%LISTOBJECT)) {
	 print STDERR "Looking into $filename\n" if ($msg >= 5);
    open(INPUT,"<", $filename) or print STDERR "ERROR: Can't open file '".$filename."\n"; #  if ($msg >= 1 )
    my $file = $LISTOBJECT{$filename};
    my $line_number = 0;

    while (<INPUT>) {
        if ($_ =~ /^[@]*[\s]*#[\s]*include[\s]*[<'"\s]([\w.\/\.]+)[>"'\s][\s]*/) {
				next if (process_file_for_include($file,$1));
        }
        next if ( $file->{EXTENSION} =~ /(c|cc|CC)$/);
        
        # FORTRAN include statement : include "..."    include ',,,"
        if ($_ =~ /^[@]*[\s]*include[\s]*[<'"\s]([\w.\/\.]+)[>"'\s][\s]*/i) {
				next if (process_file_for_include($file,$1));
        }
        # FORTRAN use statement : use yyy 
        if ($_ =~ /^[@]*[\s]*\buse[\s]+([a-z][\w]*)(,|\t| |$)/i) {
            my $modname = $1 ; $modname =~ tr/A-Z/a-z/ ; # modules names are case insensitive

            # If the module can be found, add the file to dependencies
            if (my $include_filename = search_module(\%LISTOBJECT, $modname)) {
					 ${$file->{DEPENDENCIES}}{$include_filename} = $LISTOBJECT{$include_filename} if (!exists ${$file->{DEPENDENCIES}}{$include_filename} );
					 #print STDERR "$filename +: $modname \n";
            } else {
					 push @{$file->{UNSOLVED_MODULE}}, $modname; 
					 #print STDERR "$filename -: $modname \n";
            }

        } elsif ($_ =~ /^[@]*[\s]*\buse[\s]+([a-z][\w]*)/i) {
            ${$file->{UNKOWN_USE}}{$line_number} = $_;
				#print STDERR "$filename ? \n";
        }

        # FORTRAN module declaration : module xxx
        if ($_ =~ /^[@]*[\s]*\bmodule[\s]+([a-z][\w]*)(,|\t| |$)/i) {
            my $modname = $1 ; $modname =~ tr/A-Z/a-z/ ; # modules names are case insensitive
            my $search_filename = '';

            next if $modname eq "procedure";

            #print "SEARCHING MODULE: ";
            #print "$modname\n";

            # Verifier que le nom du module n'existe pas dans un autre fichier
            if ($search_filename = search_module(\%LISTOBJECT, $modname)) { 
                print STDERR "Module ".$modname." (".$filename.") already defined in ".$search_filename."\n"; 
                next; 
            }

            # Ajouter le module dans la liste des modules associer au fichier.
            push @{$file->{ MODULE_LIST }}, $modname;

            # Recherche tous les fichiers analyser precedemment qui avait besoin de ce module la
            while(my $key = search_unsolved_module(\%LISTOBJECT, $modname)) {
                #print STDERR "unsolved module: $key".${$LISTOBJECT{$key}->{ DEPENDENCIES }}{$filename}." : $modname \n";
                # Ajouter a la liste des dependence, le fichier en cours
                ${$LISTOBJECT{$key}->{ DEPENDENCIES }}{$filename} = $file if (!exists ${$LISTOBJECT{$key}->{ DEPENDENCIES }}{$filename} );

                # Enlever le module de la liste des unsolved modules 
                $LISTOBJECT{$key}->remove_unsolved_module($modname);
            }
        } elsif ($_ =~ /^[@]*[\s]*\bmodule[\s]+/i) {
            ${$file->{ UNKOWN_MODULE }}{$line_number} = $_;
            #print STDERR "Unknown module statement: $filename: $_\n";
        }
        $line_number++;
    }
    
    $file->set_status();

    close INPUT;

}

print STDERR "Checking for Circular dependencies\n" if ($msg >= 5);
#while(my($filename, $file) = each(%LISTOBJECT)) {
for my $filename (keys %LISTOBJECT) {
	 my $file = $LISTOBJECT{$filename};
	 #print STDERR "Checking for Circular dependencies: $filename\n" if ($msg >= 5);
    my $result;
    if ($result = $file->find_depedencies($filename)) { 
        print STDERR "ERR: Circular dependencies in $filename FAILED\n";
        exit 1; 
    }
}


#
#  lists of file types FDECKS, CDECKS, ...
#
print STDERR "Listing file types FDECKS, CDECKS, ...\n" if ($msg >= 5);
reset_all_file(\%LISTOBJECT);
for $ext (keys %SRCFile::TYPE_LIST) {
    print_header(uc $ext."DECKS", "=", "");
    for (sort keys %LISTOBJECT) {
        my $file = $LISTOBJECT{$_};
		  if (lc $file->{EXTENSION} eq lc $ext) {
				if ($flat_layout) {
					 print_item("$file->{FILENAME}.$file->{EXTENSION}");
				} else {
					 print_item("$file->{FULLPATH_SRC}$file->{FILENAME}.$file->{EXTENSION}");
				}
		  }
    }
    print STDOUT "\n";
}

#
#  OBJECTS LIST
#
print STDERR "Listing OBJECTS\n" if ($msg >= 5);
print_header("OBJECTS","=","");
for (sort keys %LISTOBJECT) {
    my $file = $LISTOBJECT{$_};
    if ($file->{TYPE} eq "COMPILABLE") {
		  if ($flat_layout) {
				print_item("$file->{FILENAME}.o");
		  } else {
				print_item("$file->{FULLPATH_SRC}$file->{FILENAME}.o");
		  }
	 }
}
print STDOUT "\n";

#
#   Build dependencie rules
#
#TODO: Dependencies to Modules should be on .mod:.o not direcly on .o (.mod could have been erased accidentaly)
print STDERR "Printing dependencie rules\n" if ($msg >= 5);
for my $filename (sort keys %LISTOBJECT) {
    my $file = $LISTOBJECT{$filename};
    @current_dependencies_list = ();
    if ($file->{TYPE} eq "COMPILABLE") {
		  if ($flat_layout) {
				print_header("$file->{FILENAME}.o",":","$file->{FILENAME}.$file->{EXTENSION}");
		  } else {
				print_header("$file->{FULLPATH_SRC}$file->{FILENAME}.o",":","$filename");
		  }
        rec_print_dependencies(\%LISTOBJECT, $filename);
        print STDOUT "\n";
		  print_header("$file->{FILENAME}.o",":","$file->{FULLPATH_SRC}$file->{FILENAME}.o") if ($short_target_names and $file->{FULLPATH_SRC} and !$flat_layout);
        print STDOUT "\n";
    }
}


#
#    Print the missing module(s) / file(s) from the current tree
#
print STDERR "Includes missing from the current tree: ".join(" ",@outside_tree_list)."\n" if ($#outside_tree_list );
#TODO: do as module below, print first filename for each missing inc

%module_missing = ();
#while(my($filename, $file) = each(%LISTOBJECT)) {
for my $filename (keys %LISTOBJECT) {
    my $file = $LISTOBJECT{$filename};
    for my $module (@{$file->{UNSOLVED_MODULE}}) {
        next if ($module eq "");
		  next if (exists($module_missing_ignored{$module}));
        $module_missing{$module} = $filename if (!exists $module_missing{$module});
    }
}
if (keys %module_missing) {
	 print STDERR "Modules missing from the current tree: ";
	 while(my($module,$filename) = each(%module_missing)) {
		  print STDERR "$module ($filename) ";
	 }
	 print STDERR "\n";
}


#
#   Unknown module and use statement 
#
my $module_unknown = "";
my $use_unknown = "";
#while(my($filename, $file) = each(%LISTOBJECT)) {
for my $filename (keys %LISTOBJECT) {
    my $file = $LISTOBJECT{$filename};
    while(my($line_number,$text_line) = each(%{$file->{UNKNOWN_MODULE}})) {
		  $module_unknown .= "\t($filename) $line_number: $text_line\n";
    }
    while(my($line_number,$text_line) = each(%{$file->{UNKOWN_USE}})) {
		  $use_unknown .= "\t($filename) $line_number: $text_line\n";
    }
}
print STDERR "Unknown module statement: \n".$module_unknown if ($module_unknown);
print STDERR "Unknown use statement: \n".$use_unknown if ($use_unknown);


#
#   Export a list of produced files (.o and .mod)
#
if ($export_list) {
    open(my $EXPOUT,'>',$export_list);
    my @list_of_modules = ();
    for (keys %LISTOBJECT) {
        my $file = $LISTOBJECT{$_};
		  if ($file->{TYPE} eq "COMPILABLE") {
				if ($flat_layout) {
					 print $EXPOUT "$file->{FILENAME}.o\n";
				} else {
					 print $EXPOUT "$file->{FULLPATH_SRC}$file->{FILENAME}.o\n";
				}
		  }
        for (@{$file->{MODULE_LIST}}) {
            push @list_of_modules, $_ if $_ ne "";
        }
    }
    for (sort @list_of_modules) {
        print $EXPOUT "$_.mod\n";
    }
    close($EXPOUT);
}

## output list of produced files (.mod .o) in specified file. Else do nothing! 

# DEBUG tree
# while(my($filename, $file) = each(%LISTOBJECT) )
# {
#     $file->displ();
# }

