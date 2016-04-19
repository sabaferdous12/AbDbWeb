package antibody;

use strict;
use warnings;
use Data::Dumper;
use lib ("/home/bsm/ucbterd/perl5/lib/perl5");
use Exporter qw (import);
use List::MoreUtils qw (uniq);
use Sort::Naturally qw (nsort);

our @EXPORT_OK = qw(
get_cluster_in_array 
get_multiple_matching_files_from_directory 
print_html_header 
print_error
print_error_kabat 
print_table 
print_single_antibody 
print_multiple_antibody 
print_two_similar_antibody 
print_two_different_antibody 
check_complex 
get_resolution 
get_pdb_path
print_number
print_number_two
displayError
getNumberedFileArray
printTableKeyword
displayMessage
displayMessageSpecOrg
displayAntigenError
displayMissingFilesError);

## Reads the Redundant cluster into an array as a separate element
sub get_cluster_in_array
{
    my ($pdb_id, $filehandle) = @_;
    my @Cluster_in_array = <$filehandle>;
        
    my @clus_elem;
    my @multiple_clus;
    my @clus = grep /$pdb_id/, @Cluster_in_array;
    
    
    foreach my $elem (@clus)
    {
        chomp $elem;
        @clus_elem = split ' ', $elem;
#        my @sorted_clus = sort { substr ($a, 5) <=> substr ($b, 5) } @clus_elem;
        push (@multiple_clus, [@clus_elem]);
    }

    return @multiple_clus;
    
}

## Finds matching pdb files (to redundant cluster) from directory
sub get_multiple_matching_files_from_directory
{
    my ($dir_files_ref, $cluster_elem_ref) = @_;
    my @file_data = ();

    my @dir_files = @{$dir_files_ref};
    my @matching;
#    print "RAW CLUSTER:~~~ @$cluster_elem_ref ~~~~";
    my @sorted = sort { substr ($a, 5) <=> substr ($b, 5) } @$cluster_elem_ref;
    my @sortedClus = sort {substr ($a, 0, 5) cmp substr($b, 0, 5) } @sorted;
    
 #   print "SORTTTTTTT~~~~~~: @sort2";
    
  #  exit;
    

    foreach my $file (@sortedClus)
        {
            chomp $file;
            $file =~ tr/,//d;
            if ( $file) {
                  # To find matching element of cluster from Directory
                my @matching = grep (/^$file(.*)?\.pdb$/, @dir_files);
                push (@file_data, @matching);
                
            }
            else {
                next;               
            }
        }
    my @file_dataP = uniq (@file_data);
    # To sort numerically 4FQR_1
    my @sorted_matching = sort { substr ($a, 5) <=> substr ($b, 5) } @file_dataP;
    # To sort alphabetically 
    my @sort2 = sort {substr ($a, 0, 5) cmp substr($b, 0, 5) } @sorted_matching;
    return @sort2;
}

## Prints Header (main layout) for Html page for displaying results
sub print_html_header
{
print <<__EOF;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<link rel="stylesheet" href="http://www.bioinf.org.uk/bo.css" type="text/css" />
<title>Results</title>
</head>
<body>
<h1><b> Antibody Database </b></h1>
<h2>Results</h2>
</body>
</html>
__EOF
return 1;
}

## If antibody is not found in Database then displays this error message on new Html page
sub print_error

{
    my ($pdb_id) = @_;
print <<__EOF;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<title>Error</title>
<link rel="stylesheet" type="text/css" href="http://www.bioinf.org.uk/bo.css" />
</head>
<body>
<h1><b>Antibody Database</b></h1>
<h2>Not Found</h2>
<p>The PDB <b>$pdb_id</b> in not found in our database. 
This may be due to one of the following reasons:\
<ul>
<li>This might not be a valid PDB code for antibody structure.</li>  
<li>The PDB <b>$pdb_id</b> might have failed to process.</li>
</ul>
</p>
</body>
</html>
__EOF
    return 1;
}

## Prints table containg results (Three numbering schemes)
sub print_table
{
    
    my ($array_size, $pdb_id, $kabat_files, $chothia_files, $martin_files,
	$cgi, $less_one) = @_;
    
    my ($temp_pdb, $ext);
    print "<table class=mytable>";
    print "<tr>\n                                                             
                           
<th>Kabat Numbered</th>   
<th>Chothia Numbered</th>
<th>Martin Numbered</th>
<th>Status</th>
<th>Resolution/R-Factor</th>

</tr>";
    my ($antibody_status, $pdb_file, $resolution);
    my $pdb_bold = 0;

    for (my $i = 0; $i <= $array_size-1; $i++)
    {
        $pdb_file = substr(${$kabat_files}[$i], 0, 4);
	        
	$antibody_status = check_complex($pdb_file);
	$resolution = get_resolution($pdb_file);
               
	if (index(${$kabat_files}[$i], $pdb_id) != -1)
	{ 
## To print the Query PDB as BOLD in results table
	    
	    print "<tr>\n
 <td><a 
href=http://www.bioinf.org.uk/abs/abdb/Data/ALL_Kabat/${$kabat_files}[$i]><b>${$kabat_files}[$i]</b></a></td>
 <td><a
 href=http://www.bioinf.org.uk/abs/abdb/Data/ALL_Chothia/${$chothia_files}[$i]><b>${$chothia_files}[$i]</b></a></td>
 <td><a
 href=http://www.bioinf.org.uk/abs/abdb/Data/ALL_Martin/${$martin_files}[$i]><b>${$martin_files}[$i]</b></a></td>
<td>$antibody_status</td>
<td>$resolution</td>
</tr>";
	    
($temp_pdb, $ext)  = split ('\.', ${$kabat_files}[$i] );
$pdb_bold++;
	}
else
{ ## For printing rest of redundant antibodies in PDB
    print "<tr>\n
 <td><a href=http://www.bioinf.org.uk/abs/abdb/Data/ALL_Kabat/${$kabat_files}[$i]>${$kabat_files}[$i]</a></td>
 <td><a href=http://www.bioinf.org.uk/abs/abdb/Data/ALL_Chothia/${$chothia_files}[$i]>${$chothia_files}[$i]</a></td>
 <td><a href=http://www.bioinf.org.uk/abs/abdb/Data/ALL_Martin/${$martin_files}[$i]>${$martin_files}[$i]</a></td>
<td>$antibody_status</td>
<td>$resolution</td>
 </tr>";
}
    
    }
my $flag = 0;
if ($array_size > 2)
{
    #print_number($cgi, $temp_pdb, $pdb_id, $less_one);
    if ( ($pdb_bold == 1 )and ($array_size == 2) )
    {
	print $cgi->p("There are $less_one other antibodies that are redundant to                                 
        <b>$temp_pdb</b>.");
    }
    
    else
    {
	print $cgi->p("There are $less_one other antibodies that are redundant to                                 
        <b>$temp_pdb</b>.");
    }


}
elsif ($array_size == 2)
{
    print $cgi->p("There is $less_one other antibody that is redundant to                                       
        <b>$temp_pdb</b>.");
}

  return 1;

}

### The print sub-routines below differes only in context of displaying 
# messages for multiple or single antibody ###

## If there is only one antibody in a PDB then it just prints this message
sub print_single_antibody
{
    my ($cgi, $pdb_id) = @_;
    print $cgi->p("The antibody in PDB file <b>$pdb_id</b>,
 numbered according to 3 different numbering schemes.");
    return 1;
}

## If there are more than 2 (or redundant antibodies) antibodies in one PDB
sub print_multiple_antibody{

    my ($cgi, $pdb_id) = @_;
    print $cgi->p("The antibodies in PDB file <b>$pdb_id</b>,
 numbered according to 3 different numbering schemes.\n");
    print $cgi->p("The following table (tables) contains antibodies found "
		  ."in PDB file <b>$pdb_id</b> (in Bold), along with any ".
		  "redundant antibodies found within the database. ".
		  "Note: <b>Redundancy</b> between antibodies refers to "
		  ."100% sequence identity between heavy and light chains. "
		  ."Sequences used to test for redundancy have been extracted"
		  ." from structural data - therefore, any residues missing "
		  ."from structural data will affect sequence comparison and"
		  ." redundancy checking.");
    return 1;

}

sub print_number
{
    my ($cgi, $temp_pdb, $pdb_id, $less_one) = @_;
    print $cgi->p("There are $less_one other antibodies that are redundant to
 <b>$temp_pdb</b>.");

    return 1;
}

sub print_number_two
{
    my ($cgi, $temp_pdb, $pdb_id, $less_one) = @_;
   
 print $cgi->p("There is $less_one other antibody that is redundant to 
<b>$temp_pdb</b>.");

    return 1;

}

## If thre are 2 antibodies in a structure and there are redundant 
# e.g: 1AFV_1 and 1AFV_2
sub print_two_similar_antibody
{
    my ($cgi, $pdb_id) = @_;
    print $cgi->p("The antibodies in PDB file <b>$pdb_id</b>, numbered "
		  ."according to 3 different numbering schemes.");
    print $cgi->p("The following table contains redundant antibodies found ".
		  "in PDB file <b>$pdb_id</b>. Note: <b>Redundancy</b> "
		  ."between antibodies refers to 100% sequence identity "
		  ."between heavy and light chains. Sequences used to test ".
		  "for redundancy have been extracted from structural data "
		  ."- therefore, any residues missing from structural data "
		  ."will affect sequence comparison and redundancy checking."
	);
 
    return 1;

}

## If thre are 2 antibodies (different) for a PDB query and there are 
# redundant e.g: 1AJ7_1, 2RCS_1  
sub print_two_different_antibody
{
    my ($cgi, $pdb_id) = @_;

    print $cgi->p("The antibodies in PDB file <b>$pdb_id</b>, numbered ".
		  "according to 3 different numbering schemes.");
        print $cgi->p("The following table contains antibody found redundant "
		      ."to <b>$pdb_id</b> (in Bold). Note: <b>Redundancy</b>"
		      ." between antibodies refers to 100% sequence identity "
		      ."between heavy and light chains. Sequences used to "
		      ."test for redundancy have been extracted from "
		      ."structural data - therefore, any residues missing "
		      ."from structural data will affect sequence comparison "
		      ."and redundancy checking."
	    );

    return 1;

}

sub check_complex
{
    my ($file_name) = @_;
    my $dir = "/acrm/www/html/abs/abdb/Data";
        
    open(my $proAB,'<', "$dir/Redundant_files/Redundant_LH_Protein_Martin.txt")
	or die "Can not open file ...";
    open (my $nproAB, '<', "$dir/Redundant_files/Redundant_LH_NonProtein_Martin.txt")
        or die "Can not open file ...";
    open (my $freeAB, '<', "$dir/Redundant_files/Redundant_LH_Free_Martin.txt")
        or die "Can not open file ...";
        
    open (my $proLG, '<', "$dir/Redundant_files/Redundant_L_Protein_Martin.txt" )
        or die "Can not open file ...";
    open (my $nproLG, '<', "$dir/Redundant_files/Redundant_L_NonProtein_Martin.txt" )
        or die "Can not open file ...";
    open (my $freeLG, '<', "$dir/Redundant_files/Redundant_L_Free_Martin.txt" )
        or die "Can not open file ...";
        
    open (my $proHV, '<', "$dir/Redundant_files/Redundant_H_Protein_Martin.txt" )
        or die "Can not open file ...";
    open (my $nproHV, '<', "$dir/Redundant_files/Redundant_H_NonProtein_Martin.txt" )
        or die "Can not open file ...";
    open (my $freeHV, '<', "$dir/Redundant_files/Redundant_H_Free_Martin.txt")
        or die "Can not open file ...";
    #print "UUUUUUUUU: $file_name\n";    
    
    my @proABCluster = <$proAB>;
    my @nproABCluster = <$nproAB>; 
    my @freeAB = <$freeAB>;
    
    my @proLg = <$proLG>;
    my @nproLg = <$nproLG>;
    my @freeLg = <$freeLG>;
    
    my @proHv = <$proHV>;
    my @nproHv = <$nproHV>;
    my @freeHv = <$freeHV>;
    
    my $proAb = grep /$file_name/, @proABCluster;
    my $nproAb = grep /$file_name/, @nproABCluster;
    my $freeAb = grep /$file_name/,@freeAB;
    
    my $proLg = grep /$file_name/,@proLg;
    my $nproLg = grep /$file_name/,@nproLg;
    my $freeLg = grep /$file_name/,@freeLg;
    
    my $proHv = grep /$file_name/,@proHv;
    my $nproHv = grep /$file_name/,@nproHv;
    my $freeHv = grep /$file_name/,@freeHv;
    
    my $status;    
    if ($proAb)
    {
	$status = "Antibody Complexed (Protein)";
    }
    elsif ($nproAb)
        {
            $status = "Antibody Complexed (Non-protein)";
         }
    elsif ($freeAb)
        {
            $status = "Free Antibody";
         }
        
    elsif ($proLg ) {
        $status = "Bence-Jones Complexed (Protein)";
    }

    elsif ($nproLg ) {
        $status = "Bence-Jones Complexed (Non-Protein)";
    }

    elsif ($freeLg ) {
        $status = "Free Bence-Jones";
    }
    
    elsif ( $proHv) {
        $status = "Camelids Complexed (Protein)";
    }
    elsif ( $nproHv) {
        $status = "Camelids Complexed (Non-Protein)";
    }
    elsif ( $freeHv) {
        $status = "Free Camelids";
    }
        
      
    return $status;

}

sub get_resolution
{
    my ($pdb_file) = @_;
    my $file_path = get_pdb_path($pdb_file);
    my $prog_path = "/home/bsm/martin/bin/getresol";
    my $pdb_resol = `$prog_path $file_path`;
    my ($method, $resol) = split (",", $pdb_resol);
    
return $resol;

}


sub get_pdb_path
{
    my ($pdb_id) = @_;
    chomp $pdb_id;
    my $file_path;

    if($pdb_id=~m/[0-9A-Z]{4}/)
    {
	my $dir = '/acrm/data/pdb/';
#        my $dir = '/acrm/nouf2data1/acrm8/data/pdb/';
	my $ext = '.ent';
	my $temp = 'pdb';
	$file_path = $dir.$temp.lc($pdb_id).$ext;
    }
    else
    {
	print "This is not a valid PDB Code\n";
	exit;
    }
    
    return $file_path;
    
}

sub print_error_kabat
{
    my ($pdb_id) = @_;
    print <<__EOF;
    <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
    <html>
    <head>
    <title>Error</title>
    <link rel="stylesheet" type="text/css" href="http://www.bioinf.org.uk/bo.css" />
    </head>
    <body>
    <h1><b>Antibody Database</b></h1>
    <h2>Not Found</h2>
    <p>The antibody numbering program (Kabatnum) fails on PDB <b>$pdb_id</b>.
    <br>    
    Alternatively, <b>$pdb_id</b> can be downloaded from <a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=$pdb_id">Protein Data Bank</a>
    <br>
    <br>
    To find the complete list of antibodies with failed numbering, 
    <a href="http://www.bioinf.org.uk/abs/abdb/Data/Kabat_Failed.list">
    Click Here</a>
    </p>
    </body>
    </html>

__EOF
return 1;

}


sub displayError
{
    my ($errorStr, $errorStr2) = @_; 
print <<HTML;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<title>Error</title>
<link rel="stylesheet" type="text/css" href="http://www.bioinf.org.uk/bo.css" />
</head>
<body>
<h1>Error</h1> <p>$errorStr</p>
<p>$errorStr2</p>
</body>
</html>
HTML
exit

}

sub getNumberedFileArray
{
    my (@args) = @_; 
    my @filesArray = ();
    my @files ;
    foreach my $argPair (@args)
    {
	my @files = get_multiple_matching_files_from_directory
	    ( $argPair->[0], $argPair->[1] );
	push(@filesArray, \@files);
    }

    
    my (@kabat_files, @chothia_files, @martin_files);

    push (@kabat_files, @{$filesArray[0]});
   
    push (@chothia_files,@{$filesArray[1]});
       
    push (@martin_files, @{$filesArray[2]});
   
    return (\@kabat_files, \@chothia_files, \@martin_files);
    
}

sub printTableKeyword
    {
        my ($kabat_files, $chothia_files, $martin_files) = @_;
        my $arrSize = scalar @{$kabat_files};
        print "<table class=mytable>";
        print "<tr>\n
<th>Kabat Numbered</th>
<th>Chothia Numbered</th>
<th>Martin Numbered</th>
<th>Status</th>
<th>Resolution/R-Factor</th>
</tr>";
                 
        my ($antibody_status, $pdb_file, $resolution);
        
        for ( my $i = 0; $i <= $arrSize-1; $i++)
            {
                $pdb_file = substr(${$kabat_files}[$i], 0, 4);
                $antibody_status = check_complex($pdb_file);
                $resolution = get_resolution($pdb_file);
                
                print "<tr>\n
<td><a href=http://www.bioinf.org.uk/abs/abdb/Data/ALL_Kabat/${$kabat_files}[$i]>${$kabat_files}[$i]</a></td>
 <td><a href=http://www.bioinf.org.uk/abs/abdb/Data/ALL_Chothia/${$chothia_files}[$i]>${$chothia_files}[$i]</a></td>
 <td><a href=http://www.bioinf.org.uk/abs/abdb/Data/ALL_Martin/${$martin_files}[$i]>${$martin_files}[$i]</a></td>
<td>$antibody_status</td>
<td>$resolution</td>
 </tr>";
            }
         
    }
sub displayMessage
        {
            my ($cgi, $PDBnum, $keyword) = @_;
            # To convert the keyword into title case
            $keyword =~  s/([^\s\w]*)(\S+)/$1\u\L$2/g;
            if ( !$PDBnum) {
                print $cgi->p("There is no PDB structure found in this database for <b>$keyword</b>. This may be due to one of the following reasons:");
                print $cgi->ul(
                    $cgi->li("The structures for <b>$keyword</b> might not exists."), $cgi->li("The PDB files for <b>$keyword</b> might not have processed.") );
                exit;     
            }

            else {
                print $cgi->p("There are $PDBnum PDB structures found in this database for <b>$keyword</b>.");
            }
            
            
        }
        

sub displayMessageSpecOrg
            {
                my ($cgi, $PDBnum, $organism, $species) = @_;
                $organism =~  s/([^\s\w]*)(\S+)/$1\u\L$2/g;
                $species =~  s/([^\s\w]*)(\S+)/$1\u\L$2/g;
                 
                if  ( !$PDBnum) {
                    print $cgi->p("There is no PDB structure found in this database for antibody species <b>$organism</b> bound with antigen species <b>$species</b>. This may be due to one of the following reasons:");
                    print $cgi->ul(
                    $cgi->li("The structures for <b>$organism</b> or <b>$species</b> might not exists."), $cgi->li("The PDB files for <b>$organism</b> or  <b>$species</b> might not have processed."));
                    
                }
                else {
                    
                    print $cgi->p("There are $PDBnum PDB structures found in this database for antibody species <b>$organism</b> bound with antigen species <b>$species</b>.");
                    
                }
                
            }
            

sub displayAntigenError
{
    my ($cgi) = @_;
    print $cgi->p("<b>Note: </b>Any free antibodies appearing in the search results are due to non-antigen proteins bound to antibody (other than CDR regions). Such proteins have been removed from antibody during processing of original PDB file.");
    
}


sub displayMissingFilesError
{
    my ($cgi) = @_;
    print $cgi->p("<b>Note: </b>The missing files in following table are due to failure of one of the numbering method.");
    
}
    
1;


