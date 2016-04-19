for(my $i = 0 ; $i < @kabat ; ++$i) {
    my $kabatClusterRef = $kabat[$i];
    my $chothiaClusterRef = $chothia[$i];
    my $martinClusterRef = $martin[$i];

    my @args = ([\@dir_files_k, $kabatClusterRef], [\@dir_files_m, $martinClusterRef],
		[\@dir_files_c, $chothiaClusterRef]);

    my @filesArray = ();

    foreach my $argPair (@args) {
	my @files = get_multiple_matching_files_from_directory($argPair->[0], $argPair->[1]);
	push(@filesArray, \@files);
    }
    
    print_table($array_size, $pdb_id, @filesArray);
}
