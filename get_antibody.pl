#! /usr/bin/perl 
use warnings;
use strict;
use lib ("/home/bsm/ucbterd/perl5/lib/perl5");
use CGI;
use Data::Dumper;
use List::MoreUtils qw (uniq);
use Scalar::Util qw(looks_like_number);

use antibody qw 
    (get_cluster_in_array
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
     print_number
     displayError
     getNumberedFileArray
     get_resolution
     printTableKeyword
     displayMessage
     displayMessageSpecOrg
     displayAntigenError
     displayMissingFilesError);
my $cgi = new CGI;
print $cgi->header();

## Directory <Data> contains all the sub-directories for all the complexes
# numbered according to diffrent numbering schemes 
my $dir = "/acrm/www/html/abs/abdb/Data";
my $logDir = "/acrm/www/html/abs/abdb/Data/Martin_logs";
my $headerFile = "header.dat";
my $headerProFile = "headerProcessed2.dat";
my $antigenChainsFile = "AntigenChains.dat";
my $antibodyChainsFile = "AntibodyChains.dat";

my ($file1, $file2, $file3);
my (@kabat_files, @chothia_files, @martin_files);
my (@kabat_clus, @chothia_clus, @martin_clus);

# Reading Kabat data
open(my $CLUSTER_K,'<',
     "$dir/Redundant_files/Redundant_ALL_Kabat.txt");
opendir (my $DIR_K, "$dir/ALL_Kabat") or die
    "Can not open $dir/ALL_Kabat";
my @dir_files_k = readdir ($DIR_K);

# Reading Chothia data
open(my $CLUSTER_C,
     '<', "$dir/Redundant_files/Redundant_ALL_Chothia.txt");
opendir (my $DIR_C, "$dir/ALL_Chothia") or die
    "Can not open $dir/ALL_Chothia";
my @dir_files_c = readdir ($DIR_C);

# Reading Martin data
open(my $CLUSTER_M,'<',
     "$dir/Redundant_files/Redundant_ALL_Martin.txt");
opendir (my $DIR_M, "$dir/ALL_Martin") or die
    "Can not open $dir/ALL_Martin";
my @dir_files_m = readdir ($DIR_M);
 

my ($pdb_id, $nameAg, $speciesAg);
my ($nameAb, $speciesAb);

$pdb_id = $cgi->param("pdbid");
$pdb_id = uc($pdb_id);
$nameAg = $cgi->param("organismAg");
$nameAg = uc($nameAg );
$nameAb = $cgi->param("organismAb");
$nameAb = uc($nameAb);
$speciesAg = $cgi->param("speciesAg");
$speciesAg = uc($speciesAg);
$speciesAb = $cgi->param("speciesAb");
$speciesAb = uc($speciesAb);

# Searcg for PDB code only
if ( ($pdb_id) and ($nameAg eq "" ) and ($speciesAg eq "") )
{
    $nameAg = 0;
    $speciesAg = 0;
 
    if ($pdb_id =~ m/^[0-9]{1}[0-9A-Z]{3}$/)
    {
	@kabat_clus = get_cluster_in_array($pdb_id, $CLUSTER_K);
	@chothia_clus = get_cluster_in_array($pdb_id, $CLUSTER_C);
	@martin_clus = get_cluster_in_array($pdb_id, $CLUSTER_M);
        

        my $clusNum = @kabat_clus;
	
	if (!$clusNum)                                                         
	{                                                                         
	    open(my $KABAT_FAILED,'<', "$logDir/Kabat_Failed.list");
	    my @kabat_error = <$KABAT_FAILED>;
	    
	    if (grep (/$pdb_id/, @kabat_error) )
	    {
		print_error_kabat($pdb_id);    
		exit;
	    }
	    else 
	    {
		print_error($pdb_id);                                           
		exit;
	    }                                                                     
	}
	else
	{	
	    # If there is only one antibody in a PDB then it just prints this message 
	    print_html_header();
	    my $flag = 0;
	    my $tableNo = 0;
	    for(my $i = 0 ; $i < @kabat_clus ; ++$i) 
	    {
		my $kabatClusterRef = $kabat_clus[$i];
		my $chothiaClusterRef = $chothia_clus[$i];
		my $martinClusterRef = $martin_clus[$i];
                
		my @args = ( [\@dir_files_k, $kabatClusterRef], 
			     [\@dir_files_c, $chothiaClusterRef],
			     [\@dir_files_m, $martinClusterRef]
		    );

		my ($kabat_filesRef, $chothia_filesRef, $martin_filesRef) = 
		    getNumberedFileArray (@args);		
		
		# De-referening 
		@kabat_files = @{$kabat_filesRef};
		@chothia_files = @{$chothia_filesRef};
		@martin_files = @{$martin_filesRef};
                
                my $cluster_size = @{$kabatClusterRef};	
	
		my $less_one = $cluster_size-1;   
		
		if ($clusNum == 1)     
		{                      
		    if ($cluster_size == 1)                                           
		    {  
			print_single_antibody($cgi, $pdb_id);                         
			print_table                                            
			    ( $cluster_size, $pdb_id, \@kabat_files, \@chothia_files, 
			      \@martin_files, $cgi, $less_one );                      
		    }                 
		    
		    elsif ( ($cluster_size == 2) )                                    
		    { 
# If thre are 2 antibodies in a structure and they are redundant              
# e.g: 1AFV_1 and 1AFV_2                                                      
			if ( (index($$kabatClusterRef[0], $pdb_id) != -1) and         
			     (index($$kabatClusterRef[1], $pdb_id) != -1) )          
			{                                                             
			    print_two_similar_antibody($cgi, $pdb_id);
			    print_table                                               
				( $cluster_size, $pdb_id, 
				  \@kabat_files, \@chothia_files,\@martin_files,
				  $cgi, $less_one);                                   
			} 
# If thre are 2 antibodies (different) for a PDB query and there are redundant
# e.g: 1AJ7_1, 2RCS_1                                                         
			else                                                          
			{                                                             
			    if (!$flag)                                               
			    {                                                         
				print_two_different_antibody($cgi, $pdb_id);         
				$flag = 1;                                           
			    }                                                      
			    print_table                                            
				( $cluster_size, $pdb_id, \@kabat_files, \@chothia_files,
                              \@martin_files, $cgi, $less_one);  
			}         
		    }
		    else 
		    {
			if (!$flag)
			{
			    print_multiple_antibody($cgi, $pdb_id);
			    $flag = 1;
			}
			print_table
			    ($cluster_size, $pdb_id, \@kabat_files, \@chothia_files,
                              \@martin_files, $cgi, $less_one);
		    }
		    
		}
		elsif ($clusNum == 2)
		{
		    if ($cluster_size == 1)
		    {
			if (!$flag)      
			{ 
			    print $cgi->p("The antibodies in PDB file <b>$pdb_id</b>, numbered according to 3 different numbering schemes. Non-redundant antibodies in <b>$pdb_id</b> files are shown in separate tables.");               
			    $flag = 1;  
			} 
			print_table
			    ( $cluster_size, $pdb_id, \@kabat_files, \@chothia_files,
                              \@martin_files, $cgi, $less_one ); 
		    }       
		    elsif ($cluster_size == 2)
		    {		
# If thre are 2 antibodies in a structure and they are redundant              
# e.g: 1AFV_1 and 1AFV_2                                                      
			if ( (index($$kabatClusterRef[0], $pdb_id) != -1) and  
			     (index($$kabatClusterRef[1], $pdb_id) != -1) )       
			{                                                        
			    if (!$flag)
			    {		    
				print_two_similar_antibody($cgi, $pdb_id);  
				$flag = 1;
			    }
			    print_table                                       
				( $cluster_size, $pdb_id, \@kabat_files, \@chothia_files,
				  \@martin_files, $cgi,$less_one);
			}
			
# If thre are 2 antibodies (different) for a PDB query and there are redundant
# e.g: 1AJ7_1, 2RCS_1                                                         
			else    
			{
			    if (!$flag)  
			    {
				print_two_similar_antibody($cgi, $pdb_id);
				$flag = 1;
			    }                 	
			    print_table 
				( $cluster_size, $pdb_id, \@kabat_files, \@chothia_files,
                              \@martin_files, $cgi, $less_one);
			}
		    }
		    else
		    {
			if (!$flag)
			{
			    print_multiple_antibody($cgi, $pdb_id);
			    $flag = 1;
			}
			
			print_table
			    ($cluster_size, $pdb_id,\@kabat_files, \@chothia_files,
                              \@martin_files, $cgi, $less_one);
			
		    }
		}
		else
		{
# If there are more than 2 (or redundant antibodies) antibodies in one PDB  
		    if (!$flag)
		    {
			print_multiple_antibody($cgi, $pdb_id);
			$flag = 1;
		    }
		    
		    print_table($cluster_size, $pdb_id, \@kabat_files, \@chothia_files,
                              \@martin_files, $cgi, $less_one);                  
		    
		}
		
	    }
            
            # To check if some of the files are missing for some reason
            if ( ( @kabat_files == @chothia_files) and ( @kabat_files == @martin_files)
                                               and (@chothia_files == @martin_files) ) {
            }
            else 
            {
                displayMissingFilesError($cgi);
            }

        }
        
       
    
    }
    else
    {
	my $errorStr = "Invalid PDB code!!!";
	my $errorStr2 = "Please enter valid PDB code"; 
	displayError($errorStr, $errorStr2 );
    }
exit;     
    
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( ($nameAb) and ($nameAg) ) {
    my $errorStr = "Please search by either antibody or antigen name\n";
    displayError ($errorStr,"");
    exit;
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Search by antigen name 
if ( ($nameAg) and ($speciesAg eq "") and ($pdb_id eq "") )
{
    $speciesAg = 0;
    $pdb_id = 0;
    my $isNum = looks_like_number($nameAg); 
    if ($isNum)
	{
	    my $errorStr = "The name keyword should be a string\n";
	    displayError($errorStr, "" );
	    exit; 
	}
    else 
    {
	print_html_header();
	my @result = `grep -E MOLECULE $logDir/$headerFile | grep -E "$nameAg"`;
        my @AgMatches;
        foreach  my $org (@result) {
            my $PDB = substr ($org, 0, 4);
            my @AgChains = split (":",`grep $PDB $logDir/$antigenChainsFile`);
            my @chains = split (/\,|""/,$AgChains[1] );
                           
            foreach my $chain ( @chains) {
                chomp $chain;
                my $matchedAg =
                    `grep "^$PDB : MOLECULE : $chain :.*$nameAg" $logDir/$headerFile`;
            push (@AgMatches, $matchedAg);    
            }
    
        }
         
	my @PDBcodes = uniq (map {substr ($_, 0, 4)} @AgMatches);
        my $PDBnum;
         
        my @args = ( [\@dir_files_k, \@PDBcodes],
		     [\@dir_files_c, \@PDBcodes],
		     [\@dir_files_m, \@PDBcodes]
	    ); 
	
	my ($kabat_filesRef, $chothia_filesRef, $martin_filesRef) =
	    getNumberedFileArray (@args);
        $PDBnum = scalar @{$kabat_filesRef};
        
        displayMessage ($cgi, $PDBnum, $nameAg); 
        if ( !$PDBnum) {
            exit;

        }
        # De-referening                                                                                                          
	@kabat_files = @{$kabat_filesRef};
	@chothia_files = @{$chothia_filesRef};
	@martin_files = @{$martin_filesRef};

        printTableKeyword(\@kabat_files, \@chothia_files, \@martin_files);


          # To check if some of the files are missing for some reason
        if  ( ( @kabat_files == @chothia_files) and ( @kabat_files == @martin_files)
                  and (@chothia_files == @martin_files) ) {
        }
        else {
            displayMissingFilesError($cgi);
        }
        
    }
    
    
    exit; 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Search by antibody name (keyword)
if ( ($nameAb) and ($speciesAg eq "") and ($pdb_id eq "") )
{
    $speciesAg = 0;
    $pdb_id = 0;
    my $isNum = looks_like_number($nameAb); 
    if ($isNum)
	{
	    my $errorStr = "The name keyword should be a string\n";
	    displayError($errorStr, "" );
	    exit; 
	}
    else 
    {
	print_html_header();
	my @result = `grep -E MOLECULE $logDir/$headerFile | grep -E "$nameAb"`;
        my @AbMatches;
        foreach  my $org (@result) {
            my $PDB = substr ($org, 0, 4);
            my @AbChains = split (":",`grep $PDB $logDir/$antibodyChainsFile`);
            my @chains = split (/\,|""/,$AbChains[1] );
                           
            foreach my $chain ( @chains) {
                chomp $chain;
                my $matchedAb =
                    `grep "^$PDB : MOLECULE : $chain :.*$nameAb" $logDir/$headerFile`;
            push (@AbMatches, $matchedAb);    
            }
    
        }
         
	my @PDBcodes = uniq (map {substr ($_, 0, 4)} @AbMatches);
        my $PDBnum;
         
        my @args = ( [\@dir_files_k, \@PDBcodes],
		     [\@dir_files_c, \@PDBcodes],
		     [\@dir_files_m, \@PDBcodes]
	    ); 
	
	my ($kabat_filesRef, $chothia_filesRef, $martin_filesRef) =
	    getNumberedFileArray (@args);
        $PDBnum = scalar @{$kabat_filesRef};
        
        displayMessage ($cgi, $PDBnum, $nameAb); 
        if ( !$PDBnum) {
            exit;

        }
        # De-referening                                                                                                         
	@kabat_files = @{$kabat_filesRef};
	@chothia_files = @{$chothia_filesRef};
	@martin_files = @{$martin_filesRef};

        printTableKeyword(\@kabat_files, \@chothia_files, \@martin_files);


        # To check if some of the files are missing for some reason
        if  ( ( @kabat_files == @chothia_files) and ( @kabat_files == @martin_files)
                  and (@chothia_files == @martin_files) ) {
        }
        else {
            displayMissingFilesError($cgi);
        }
 
       
    }
    
    exit; 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Search for antibody and antigen species
if ( ( ($speciesAb) and ($speciesAg) and ($pdb_id eq "") ))
{
    $pdb_id = 0;    
    print_html_header();       
####################    

    my @PDBcodes = getTwoSpecies($speciesAb, $speciesAg, $antibodyChainsFile, $antigenChainsFile);
    
    my $PDBnumS;

    my @args = ( [\@dir_files_k, \@PDBcodes],
		 [\@dir_files_c, \@PDBcodes],
		 [\@dir_files_m, \@PDBcodes]
	);

    my ($kabat_filesRef, $chothia_filesRef, $martin_filesRef) =
	getNumberedFileArray (@args);
    $PDBnumS = scalar @{$kabat_filesRef};
    displayMessageSpecOrg ($cgi, $PDBnumS, $speciesAb, $speciesAg);
    if ( !$PDBnumS ) {
        exit;
    }
# De-referening                                                                                           
    @kabat_files = @{$kabat_filesRef};
    @chothia_files = @{$chothia_filesRef};
    @martin_files = @{$martin_filesRef};
    
    printTableKeyword(\@kabat_files, \@chothia_files, \@martin_files);
    displayAntigenError($cgi);

    # To check if some of the files are missing for some reason
    if  ( ( @kabat_files == @chothia_files) and ( @kabat_files == @martin_files)
              and (@chothia_files == @martin_files) ) {
    }
    else {
        displayMissingFilesError($cgi);
    }
    exit;     
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Search for antigen species only
if ( ($speciesAg) and ($nameAg eq "") and ($pdb_id eq "") )
{
    $nameAg = 0;
    $pdb_id = 0;    
    print_html_header();       
    
    my @PDBcodesS = getSpecies($speciesAg, $antigenChainsFile);
    
    my $PDBnumS;
   
    my @args = ( [\@dir_files_k, \@PDBcodesS],
		 [\@dir_files_c, \@PDBcodesS],
		 [\@dir_files_m, \@PDBcodesS]
	);

    my ($kabat_filesRef, $chothia_filesRef, $martin_filesRef) =
	getNumberedFileArray (@args);
    $PDBnumS = scalar @{$kabat_filesRef};
    displayMessage ($cgi, $PDBnumS, $speciesAg);
    if ( !$PDBnumS ) {
        exit;
    }
# De-referening                                                                                           
    @kabat_files = @{$kabat_filesRef};
    @chothia_files = @{$chothia_filesRef};
    @martin_files = @{$martin_filesRef};
    
    printTableKeyword(\@kabat_files, \@chothia_files, \@martin_files);
    displayAntigenError($cgi);

    # To check if some of the files are missing for some reason
    if  ( ( @kabat_files == @chothia_files) and ( @kabat_files == @martin_files)
              and (@chothia_files == @martin_files) ) {
    }
    else {
        displayMissingFilesError($cgi);
    }

    
    exit;     
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Search for antibody species 
if ( ($speciesAb) and ($nameAg eq "") and ($pdb_id eq "") )
{
    $nameAg = 0;
    $pdb_id = 0;    
    print_html_header();       
    
    my @PDBcodesS = getSpecies ($speciesAb, $antibodyChainsFile);
    
    my $PDBnumS;
   
    my @args = ( [\@dir_files_k, \@PDBcodesS],
		 [\@dir_files_c, \@PDBcodesS],
		 [\@dir_files_m, \@PDBcodesS]
	);

    my ($kabat_filesRef, $chothia_filesRef, $martin_filesRef) =
	getNumberedFileArray (@args);
    $PDBnumS = scalar @{$kabat_filesRef};
    displayMessage ($cgi, $PDBnumS, $speciesAb);
    if ( !$PDBnumS ) {
        exit;
    }
# De-referening                                                                                           
    @kabat_files = @{$kabat_filesRef};
    @chothia_files = @{$chothia_filesRef};
    @martin_files = @{$martin_filesRef};
    
    printTableKeyword(\@kabat_files, \@chothia_files, \@martin_files); 

    # To check if some of the files are missing for some reason
    if  ( ( @kabat_files == @chothia_files) and ( @kabat_files == @martin_files)
              and (@chothia_files == @martin_files) ) {
    }
    else {
        displayMissingFilesError($cgi);
    }

    exit;     
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This is non-functional piece of code.. and needs to be removed
=pod
=cut

    
if (( $pdb_id) and (($nameAg) or ($speciesAg)) )
    {
        my $errorStr = "Please search by either PDB code or Keyword";
        displayError ($errorStr,"");
    }

else
    {
        my $errorStr = "Invalid PDB code or keyword!!!";
        my $errorStr2 = "Please enter valid PDB code or keyword as string";
        displayError ($errorStr, $errorStr2);
    }

sub getSpecies
    {
        my ($species, $chainsFile) = @_;
        
        my @resultS = `grep SPECIES $logDir/$headerFile | grep "$species"`;
        my @Matches;
        
        foreach  my $org (@resultS) {
            my $PDB = substr ($org, 0, 4);
            
            my @AbOrAgChains = split (":",`grep $PDB $logDir/$chainsFile`);
            my @chains = split (/\,|""/,$AbOrAgChains[1] );
                
            foreach my $chain ( @chains) {
                chomp $chain;
                my $matchedAbOrAg =
                    `grep "^$PDB : SPECIES  : $chain : $species" $logDir/$headerFile`;
                push (@Matches, $matchedAbOrAg);
            }
 
        }
        my @PDBcodesS = uniq (map {substr ($_, 0, 4)} @Matches);

        return @PDBcodesS;
        
    }
    
sub getTwoSpecies
        {
            my ($speciesAb, $speciesAg, $AbchainsFile, $AgchainsFile) = @_;

            # Grep all lines for searched species
            my @resultSAb = `grep SPECIES $logDir/$headerFile | grep "$speciesAb"`;
            my @Matches;
            
            foreach  my $org (@resultSAb) {
                my $PDB = substr ($org, 0, 4);
                # Look for chain labels for species PDB
                my @AbChains = split (":",`grep $PDB $logDir/$AbchainsFile`);
                my @AgChains = split (":",`grep $PDB $logDir/$AgchainsFile`);

                # In case Ag has no chains, it will have space
                # To remove space
                if (  $AgChains[1] =~ m/\s+/) {
                   $AgChains[1] =~ s/\s+$//;
                 #   next;
               }

                next if (  !$AgChains[1]);
                                
                my @ABchains = split (/\,|""/,$AbChains[1] );
                my @AGchains = split (/\,|""/,$AgChains[1] );

                
                my @twoSpeciesRes;
                my (@ABSpecies, @AGSpecies);

                for (my $i = 0; $i< scalar @ABchains; $i++ )
                    {

                        @ABSpecies = `grep "^$PDB : SPECIES  : $ABchains[$i] : $speciesAb" $logDir/$headerFile`;
                        @AGSpecies = `grep "^$PDB : SPECIES  : $AGchains[$i] : $speciesAg" $logDir/$headerFile`;

                        if ( ( @ABSpecies) and (@AGSpecies) )  {
                            push (@Matches, $PDB);
                        }

                    }
                
            }

            
            my @PDBcodesS = uniq (map {substr ($_, 0, 4)} @Matches);
            
            return @PDBcodesS;
            
        }
        
