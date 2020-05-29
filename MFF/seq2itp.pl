#!/usr/bin/perl
# This script generates a MARTINI coarse-grained topology for proteins based 
# on the aminoacid sequence and secondary structure.
#
# For more information type './seq2itp.pl -h'
# 
# Written by: Senthil K. Kandasamy
# Revised by: Martti Louhivuori (m.j.louhivuori@rug.nl)
# 
# Version 1.1.2
# 15.8.2008
#   Fixed side-chain bead selection and angles for tryptophans.
#
# Version 1.1.1
# ---
# 18.6.2008
#   Fixed back-bone bead selection for prolines.
#
# Version 1.1
# ---
# 27.3.2008
#   Added support for using elastic networks instead of dihedrals for extended 
#   segments. -- ML
# 21.5.2008
#   Fixed inconsistencies with the published article. E.g. use the weakest 
#   backbone bond force constant in a segment and handle prolines in helices 
#   correctly. -- ML
# 4.6.2008
#   Changed coil force constant from 200 to 400 for improved stability. -- ML
# 16.6.2008
#   Added documentation, rewrote command-line argument handling and output, 
#   clarified code. -- ML
#
use warnings;
use Getopt::Long;

$my_name = "seq2itp.pl";
$version = "1.1.2";
$martini = "2.1";
$terminal_width = 78;


#---------------------------------------------------------------------------#
# GLOBAL HASHES
# 

# main parameters
our %parameter_hash = (
BEAD_BBN_HLX => N0,
BEAD_BBN_COI => P5,
BEAD_BBN_EXT => Nda,
BEAD_BBN_TRN => Nda,
BEAD_BBN_BET => Nda,
BEAD_BBN_BND => P5,
BEAD_BBN_5HX => N0,
BEAD_BBN_3HX => N0,
BEAD_BBN_NTR => Qd,
BEAD_BBN_CTR => Qa,
BEAD_BBH_CN1 => Nd,
BEAD_BBH_CN2 => Na,
BEAD_BBH_CN3 => Nda,
BEAD_ALA_HLX => C5,
BEAD_ALA_EXT => N0,
BEAD_ALA_COI => P4,
BEAD_ALA_TRN => N0,
BEAD_ALA_BND => P4,
BEAD_GLY_HLX => N0,
BEAD_GLY_EXT => Nda,
BEAD_GLY_COI => P5,
BEAD_GLY_TRN => Nda,
BEAD_GLY_BND => P5,
BEAD_PRO_HLX => C5,
BEAD_PRO_EXT => N0,
BEAD_PRO_COI => Na,
BEAD_PRO_TRN => N0,
BEAD_PRO_BND => Na,
BEAD_SC1_ASN => P5,
BEAD_SC1_ASP => Qa,
BEAD_SC1_GLU => Qa,
BEAD_SC1_GLN => P4,
BEAD_SC1_VAL => AC2,
BEAD_SC1_LEU => AC1,
BEAD_SC1_ILE => AC1,
BEAD_SC1_MET => C5,
BEAD_SC1_THR => P1,
BEAD_SC1_SER => P1,
BEAD_SC1_CYS => C5,
BEAD_SC1_PRO => AC2,
BEAD_SC1_LYS => C3,
BEAD_SC2_LYS => Qd,
BEAD_SC1_ARG => N0,
BEAD_SC2_ARG => Qd,
BEAD_SC1_HIS => SC4,
BEAD_SC2_HIS => SP1,
BEAD_SC3_HIS => SP1,
BEAD_SC1_TYR => SC4,
BEAD_SC2_TYR => SC4,
BEAD_SC3_TYR => SP1,
BEAD_SC1_PHE => SC4,
BEAD_SC2_PHE => SC4,
BEAD_SC3_PHE => SC4,
BEAD_SC1_TRP => SC4,
BEAD_SC2_TRP => SP1,
BEAD_SC3_TRP => SC4,
BEAD_SC4_TRP => SC4,
BNLN_BBN_HLX => 0.35,
BNLN_BBN_COI => 0.35,
BNLN_BBN_EXT => 0.35,
BNLN_BBN_TRN => 0.35,
BNLN_BBN_BET => 0.35,
BNLN_BBN_BND => 0.35,
BNLN_BBN_5HX => 0.35,
BNLN_BBN_3HX => 0.35,
BNLN_SC1_ASN => 0.32,
BNLN_SC1_ASP => 0.32,
BNLN_SC1_GLU => 0.40,
BNLN_SC1_GLN => 0.40,
BNLN_SC1_VAL => 0.265,
BNLN_SC1_LEU => 0.33,
BNLN_SC1_ILE => 0.31,
BNLN_SC1_MET => 0.40,
BNLN_SC1_THR => 0.26,
BNLN_SC1_SER => 0.25,
BNLN_SC1_CYS => 0.31,
BNLN_SC1_PRO => 0.3,
BNLN_SC1_LYS => 0.33,
BNLN_SC2_LYS => 0.28,
BNLN_SC1_ARG => 0.33,
BNLN_SC2_ARG => 0.34,
BNLN_SC1_HIS => 0.32,
BNLN_TR1_HIS => 0.27,
BNLN_TR2_HIS => 0.27,
BNLN_TR3_HIS => 0.27,
BNLN_SC1_TYR => 0.32,
BNLN_TR1_TYR => 0.27,
BNLN_TR2_TYR => 0.27,
BNLN_TR3_TYR => 0.27,
BNLN_SC1_PHE => 0.31,
BNLN_TR1_PHE => 0.27,
BNLN_TR2_PHE => 0.27,
BNLN_TR3_PHE => 0.27,
BNLN_SC1_TRP => 0.3,
BNLN_TR1_TRP => 0.27,
BNLN_TR2_TRP => 0.27,
BNLN_TR3_TRP => 0.27,
BNLN_TR4_TRP => 0.27,
BNLN_TR5_TRP => 0.27,
BNLN_BRD_CYS => 0.39,
BNLN_HBN_HLX => 0.61,
BNKB_BBN_HLX => 1250,
BNKB_BBN_COI => 400,
BNKB_BBN_EXT => 1250,
BNKB_BBN_TRN => 500,
BNKB_BBN_BET => 1250,
BNKB_BBN_BND => 400,
BNKB_BBN_5HX => 1250,
BNKB_BBN_3HX => 1250,
BNKB_SC1_ASN => 5000,
BNKB_SC1_ASP => 7500,
BNKB_SC1_GLU => 5000,
BNKB_SC1_GLN => 5000,
BNKB_SC1_VAL => 7500,
BNKB_SC1_LEU => 7500,
BNKB_SC1_ILE => 7500,
BNKB_SC1_MET => 2500,
BNKB_SC1_THR => 7500,
BNKB_SC1_SER => 7500,
BNKB_SC1_CYS => 7500,
BNKB_SC1_PRO => 7500,
BNKB_SC1_LYS => 5000,
BNKB_SC2_LYS => 5000,
BNKB_SC1_ARG => 5000,
BNKB_SC2_ARG => 5000,
BNKB_SC1_HIS => 7500,
BNKB_TR1_HIS => 1000,
BNKB_TR2_HIS => 1000,
BNKB_TR3_HIS => 1000,
BNKB_SC1_TYR => 5000,
BNKB_TR1_TYR => 1000,
BNKB_TR2_TYR => 1000,
BNKB_TR3_TYR => 1000,
BNKB_SC1_PHE => 7500,
BNKB_TR1_PHE => 1000,
BNKB_TR2_PHE => 1000,
BNKB_TR3_PHE => 1000,
BNKB_SC1_TRP => 5000,
BNKB_TR1_TRP => 1000,
BNKB_TR2_TRP => 1000,
BNKB_TR3_TRP => 1000,
BNKB_TR4_TRP => 1000,
BNKB_BRD_CYS => 5000,
BNKB_HBN_HLX => 1250,
BNKB_RNG_EQL => 1000,
ANGL_BBN_HLX => 96,
ANGL_BBN_COI => 127,
ANGL_BBN_EXT => 134,
ANGL_BBN_TRN => 100,
ANGL_BBN_BET => 134,
ANGL_BBN_BND => 130,
ANGL_BBN_5HX => 96,
ANGL_BBN_3HX => 96,
ANGL_BBN_HPR => 98,
ANGL_SC1_ASN => 100,
ANGL_SC1_ASP => 100,
ANGL_SC1_GLU => 100,
ANGL_SC1_GLN => 100,
ANGL_SC1_VAL => 100,
ANGL_SC1_LEU => 100,
ANGL_SC1_ILE => 100,
ANGL_SC1_MET => 100,
ANGL_SC1_THR => 100,
ANGL_SC1_SER => 100,
ANGL_SC1_CYS => 100,
ANGL_SC1_PRO => 100,
ANGL_SC1_LYS => 100,
ANGL_SC2_LYS => 180,
ANGL_SC1_ARG => 100,
ANGL_SC2_ARG => 180,
ANGL_SC1_HIS => 100,
ANGL_SC2_HIS => 150,
ANGL_SC3_HIS => 150,
ANGL_SC1_TYR => 100,
ANGL_SC2_TYR => 150,
ANGL_SC3_TYR => 150,
ANGL_SC1_PHE => 100,
ANGL_SC2_PHE => 150,
ANGL_SC3_PHE => 150,
ANGL_SC1_TRP => 100,
ANGL_SC2_TRP => 210,
ANGL_SC3_TRP => 90,
ANKB_BBN_HLX => 700,
ANKB_BBN_COI => 25,
ANKB_BBN_EXT => 25,
ANKB_BBN_TRN => 25,
ANKB_BBN_BET => 25,
ANKB_BBN_BND => 25,
ANKB_BBN_5HX => 700,
ANKB_BBN_3HX => 700,
ANKB_BBN_HPR => 100,
ANKB_SC1_ASN => 25,
ANKB_SC1_ASP => 25,
ANKB_SC1_GLN => 25,
ANKB_SC1_GLU => 25,
ANKB_SC1_VAL => 25,
ANKB_SC1_LEU => 25,
ANKB_SC1_ILE => 25,
ANKB_SC1_MET => 25,
ANKB_SC1_THR => 25,
ANKB_SC1_SER => 25,
ANKB_SC1_CYS => 25,
ANKB_SC1_LYS => 25,
ANKB_SC2_LYS => 25,
ANKB_SC1_PRO => 25,
ANKB_SC1_ARG => 25,
ANKB_SC2_ARG => 25,
ANKB_SC1_HIS => 25,
ANKB_SC2_HIS => 50,
ANKB_SC3_HIS => 50,
ANKB_SC1_TYR => 25,
ANKB_SC2_TYR => 50,
ANKB_SC3_TYR => 50,
ANKB_SC1_PHE => 25,
ANKB_SC2_PHE => 50,
ANKB_SC3_PHE => 50,
ANKB_SC1_TRP => 25,
ANKB_SC2_TRP => 50,
ANKB_SC3_TRP => 50,
DIAN_BBN_HLX => -120,
DIAN_BBN_EXT => 0,
DIAN_BBN_COI => 0,
DIAN_RG1_PHE => 0,
DIAN_RG1_HIS => 0,
DIAN_RG1_TYR => 0,
DIAN_RG1_TRP => 0,
DIAN_RG2_TRP => 0,
DIKB_BBN_HLX => 400,
DIKB_BBN_EXT => 10,
DIKB_BBN_COI => 10,
DIKB_RG1_PHE => 50,
DIKB_RG1_HIS => 50,
DIKB_RG1_TYR => 50,
DIKB_RG1_TRP => 50,
DIKB_RG2_TRP => 200,
DIKB_HLX_PRO => 100,
YES_BBN_CONST15 => "no",
YES_BBN_PROPDIH => "yes",
BNLN_ELASTIC_SHORT => 0.64,
BNLN_ELASTIC_LONG => 0.97,
BNKB_ELASTIC => 2500,
ITV_BONDS => "constraints",
RING_BONDS => "constraints",
Q_TERMINI => "CHARGED",
);

# aminoacids
our %aa_hash = (
	G=>"GLY", A=>"ALA", D=>"ASP", N=>"ASN", E=>"GLU", Q=>"GLN", V=>"VAL", 
	L=>"LEU", I=>"ILE", M=>"MET", T=>"THR", S=>"SER", C=>"CYS", K=>"LYS", 
	R=>"ARG", H=>"HIS", F=>"PHE", P=>"PRO", W=>"TRP", Y=>"TYR");

# secondary structure elements
our %ss_hash = (
	C=>"COI", H=>"HLX", S=>"BND", T=>"TRN", B=>"BET", G=>"3HX", I=>"5HX", 
	E=>"EXT");
our %ss_hash_inv = reverse %ss_hash;

# charges
our %q_hash = (
	G => 0.0, A => 0.0, D =>-1.0, N => 0.0, E =>-1.0, Q => 0.0, V => 0.0, 
	L => 0.0, I => 0.0, M => 0.0, T => 0.0, S => 0.0, C => 0.0, K => 1.0, 
	R => 1.0, H => 0.0, F => 0.0, P => 0.0, W => 0.0, Y => 0.0);

# number of CG beads
our %numbeads_hash = (
	G=>1, A=>1, D=>2, N=>2, E=>2, Q=>2, V=>2, L=>2, I=>2, M=>2, T=>2, S=>2, 
	C=>2, K=>3, R=>3, H=>4, F=>4, P=>2, W=>5, Y=>4);


#---------------------------------------------------------------------------#
# SUB ROUTINES
# 

#
# Echo text snippets formatted in man-style.
#   i       -- (int)    indent level
#   txt     -- (string) paragraph text to be formatted
#   postfix -- (string) end paragraph with this string
#
sub help_echo {
	local($i, $txt, $postfix) = ($_[0], $_[1], $_[2]);
	local($indent) = '';
	while (length $indent < $i * 2) {
		$indent .= '  ';
	}
	local($width) = $terminal_width - length $indent;

	local(@words) = split / /, $txt;
	local($line) = '';
	while (scalar @words) {
		local($word) = shift @words;
		if ((length $line.$word) >= $width) {
			print STDERR $indent, substr($line, 0, -1), "\n";
			$line = '';
		}
		$line .= $word.' ';
	}
	print STDERR $indent, substr($line, 0, -1), "\n$postfix";
}

#
# Echo help paragraph.
#   i   -- (int)    indent level
#   txt -- (string) text paragraph
# 
sub help {
	local($i, $txt) = ($_[0], $_[1]);
	help_echo $i, $txt, "\n";
}

#
# Echo help paragraph without an ending line-feed.
#   i   -- (int)    indent level
#   txt -- (string) text paragraph
# 
sub help_nolf {
	local($i, $txt) = ($_[0], $_[1]);
	help_echo $i, $txt, "";
}

#
# Echo man-page header, i.e. NAME and SYNOPSIS
# 
sub echo_man_head {
	help 0, "NAME";
	help_nolf 1, "$my_name (v.$version)";
	help 2, "Generate a MARTINI $martini coarse-grained topology for proteins.";
	help 0, "SYNOPSIS";
	help 1, "$my_name [options] protein.seq protein.ssd";
}

#
# Echo man-page
# 
sub echo_man {
	echo_man_head;
	help 0, "DESCRIPTION";
	help 1, "The coarse-grained MARTINI topology is generated from the aminoacid sequence and secondary structure information by selecting matching bead types, correct force constants and appropriate constraints. For further information see Monticelli et al. (2008) J Chem Theory and Comput, 4: 819-834.";
	help 1, "The bead types are selected in the following manner. For single-bead residues (i.e. ALA and GLY) and prolines the bead type depends on the secondary structure. For residues with multiple beads the backbone bead is selected according to the secondary structure whereas the side-chain bead(s) depend on the aminoacid.";
	help 1, "Force constants follow those presented in Monticelli et al. (2008) with the exception of coil angle force constant that is slightly increased to improve the stability of the model in some cases. For the same reason local elastic networks are used by default instead of dihedrals to restrain extended regions. The following special rules are also applied. For backbone bonds between two residues that have a different secondary structure, the weakest force constant is used for the bond. The backbone-backbone-backbone angles are determined by the secondary structure of the middle residue unless all three residues are either in coil, extended, turn or bend conformation in which case coil parameters are used. If the middle residue is a proline in helical conformation then a lowered force constant is used to reflect the helix-breaking propensity of prolines. The parameters of helix backbone dihedrals are set based on the first central residue unless either of the central residues is a proline in which case a lowered force constant is used once again.";
	help 0, "ARGUMENTS"; 
	help 1, "{-s/--sequence} protein.seq";
	help 2, "The aminoacid sequence of the protein, in FASTA format.";
	help 1, "{-2/--ssd/--structure} protein.ssd";
	help 2, "The secondary structure of the residues in SSDUMP format (see for example the program do_dssp in GROMACS).";
	help 0, "OPTIONS";
	help 1, "-h/-?/--help";
	help 2, "Help!";
	help 1, "-p/--par/--custom-parameters FILENAME";
	help 2, "Use custom parameters defined in FILENAME instead of the default ones. (Hint: option -out may be used to get a template for the file)";
	help 1, "-o/--output-parameters FILENAME";
	help 2, "Write all the parameters used into FILENAME.";
	help 1, "-c/--cysteines FILENAME";
	help 2, "Build disulphide bridges between cysteines defined in FILENAME.";
	help 1, "-e/--elastic";
	help 2, "Use local elastic networks to restraint extended regions. (default)";
	help 1, "-n/--noelastic/--dihedrals";
	help 2, "Use dihedrals to restraint extended regions instead of local elastic networks.";
	help 1, "-t/--itp/--topology protein.itp";
	help 2, "Write topology to this file (default: STDOUT)";
	help 1, "-v/--verbose";
	help 2, "Be loud and noisy.";
	help 1, "-d/--debug";
	help 2, "Debug.";
	help 0, "FILES";
	help 1, "Aminoacid sequence";
	help 2, "The FASTA format used for the aminoacid sequence file has two main elements: chain identifiers and aminoacid sequences. Each chain is defined as a separate entity containing one line with the chain identifier and a second line with the aminoacid sequence of the residues in the chain. E.g.";
	help_nolf 2, ">1PGA:A";
	help 2, "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
	help 2, "And for a protein that has multiple chains this becomes e.g.";
	help_nolf 2, ">MSCL:A";
	help_nolf 2, "MLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGILRIG";
	help_nolf 2, "IGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDL";
	help_nolf 2, "LAQTN";
	help_nolf 2, ">MSCL:B";
	help_nolf 2, "MLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGILRIG";
	help_nolf 2, "IGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDL";
	help_nolf 2, "LAQTN";
	help_nolf 2, ">MSCL:C";
	help_nolf 2, "MLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGILRIG";
	help_nolf 2, "IGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDL";
	help_nolf 2, "LAQTN";
	help_nolf 2, ">MSCL:D";
	help_nolf 2, "MLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGILRIG";
	help_nolf 2, "IGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDL";
	help_nolf 2, "LAQTN";
	help_nolf 2, ">MSCL:E";
	help_nolf 2, "MLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGILRIG";
	help_nolf 2, "IGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDL";
	help 2, "LAQTN";
	
	help 1, "Secondary structure";
	help 2, "SSDUMP format is used to define the secondary structure of the residues. The first line has the number of residues in the whole protein and the second line has either H, E or C for each residue marking, accordingly, either helix, extended or coil structure. E.g.";
	help_nolf 2, "56";
	help 2, "EEEEEEEEECCEEEEEEEEECCHHHHHHHHHHHHHHCCCCCEEEEECCCEEEEEEE";

	help 0, "EXAMPLES";
	help 1, "Generate a topology for protein G subunit B1 (PDB: 1pga).";
	help_nolf 2, "1. Get the sequence, e.g.";
	help_nolf 3, "> wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt";
	help 3, "> grep -A 1 1pga pdb_seqres.txt > 1pga.seq";
	help 2, "2. Create a SSDUMP file that defines the secondary structure (e.g. with the Gromacs tool do_dssp).";
	help_nolf 2, "3. Generate the topology";
	help 3, "> ./seq2itp.pl 1pga.seq 1pga.ssd > 1pga.itp";
	help 1, "If dihedrals are desired instead of local elastic networks the command would be";
	help 2, "> ./seq2itp.pl --noelastic 1pga.seq 1pga.ssd > 1pga.itp";
	help 0, "AUTHOR";
	help 1, "This script was originally written by Senthil K. Kandasamy and subsequently revised by Martti Louhivuori (m.j.louhivuori\@rug.nl).";
}

#
# Echo how-to see full documentation.
#
sub echo_help_hint {
	print STDERR "For more information, use -h / --help.\n";
}

#
# Report an error.
#   msg -- (string) error message
#
sub error {
	print STDERR "Error: $_[0]\n";
	echo_help_hint;
	exit;
}

#
# Create a single CG bead.
#   id    -- (int) atom index
#   resid -- (int) residue index
#   type  -- (str) bead type, i.e. single / backbone / SC1 etc.
#
sub atom_bead {
	my ($i, $r, $type) = @_;
	if ($DEBUG) {
		print STDERR "DEBUG: atom_bead() < $i, $r, $type\n";
	}
	my $aa = $sequence[$r];
	my $ss = $structure[$r];
	my $bead = "";
	my $name = "";
	if ($type eq 'single') {
		$bead = $parameter_hash{"BEAD_$aa_hash{$aa}_$ss_hash{$ss}"};
		$name = "B$ss" . $bead; 
	} elsif ($type eq 'proline') {
		$bead = $parameter_hash{"BEAD_PRO_$ss_hash{$ss}"};
		$name = "B$ss" . $bead;
	} elsif ($type eq 'backbone') {
		$bead = $parameter_hash{"BEAD_BBN_$ss_hash{$ss}"};
		$name = "B$ss" . $bead;
	} else {
		$bead = $parameter_hash{"BEAD_$type" . "_" . $aa_hash{$aa}};
		$name = "S$ss" . $bead;
	}
	my %atom = (
		id => $i + 1,
		bead => $bead,
		resid => $r + 1,
		resname => $aa_hash{$aa},
		name => $name,
		charge => 0.0, 
		conf => $ss_hash{$ss}
	);
	if ($DEBUG) {
		print STDERR "DEBUG atom_bead() > ";
		print STDERR map { "$_=>$atom{$_}, " } keys %atom;
		print STDERR "\n";
	}
	return %atom;
}

#
# Make atom charged.
#   i -- (int) atom index
#   
sub charge_atom {
	my $i = shift;
	$atoms[$i]{charge} = $q_hash{ $sequence[$atoms[$i]{resid} - 1] };
}

#
# Update bead name.
#   i -- (int) atom index
#   
sub update_bead_name {
	my $i = shift;
	my $ss = $ss_hash_inv{$atoms[$i]{conf}};
	my $old = $atoms[$i]{name};
	my $prefix = substr $atoms[$i]{name}, 0, 1;
	$atoms[$i]{name} = $prefix . $ss . $atoms[$i]{bead};
	if ($DEBUG) {
		print STDERR "DEBUG: update_bead_name() < $i\n";
		my $new = $atoms[$i]{name};
		print STDERR "DEBUG: updated bead name from $old to $new\n";
	}
}

#
# Create CG atoms for a residue.
#   count -- (int) number of atoms already created
#   resid -- (int) residue index
#
sub create_atoms {
	my $n = $numbeads_hash{$sequence[$_[1]]};
	# single bead?
	if ($n == 1) {
		push @atoms, { (atom_bead($_[0], $_[1], 'single')) };
	} else {
		if ( $sequence[$_[1]] eq 'P' ) {
			push @atoms, { (atom_bead($_[0], $_[1], 'proline')) };
		} else {
			push @atoms, { (atom_bead($_[0], $_[1], 'backbone')) };
		}
		my $i = 1;
		while ($i < $n) {
			push @atoms, { (atom_bead($_[0] + $i, $_[1], "SC$i")) };
			$i++;
		}
		charge_atom($#atoms);
	}
	$backbone[$_[1]] = $_[0] + 1;
	return $n;
}

#
# Echo a bond definition.
#   i      -- (int)   ID of 1st atom
#   j      -- (int)   ID of 2nd atom
#   length -- (float) bond length
#   force  -- (int)   bond force constant
#   name   -- (str)   desciption of bond
#
sub echo_bond {
	my ($i, $j, $length, $force, $name) = @_;
	$length = $parameter_hash{"BNLN_$length"};
	$force = $parameter_hash{"BNKB_$force"};
	printf("%d\t%d\t1\t%1.3f\t%d\t;\t%s\n", $i, $j, $length, $force, $name);
}

#
# Echo a sidechain bond definition.
#   i    -- (int) index of backbone atom
#   s    -- (int) ID shift of 1st atom
#   t    -- (int) ID shift of 2nd atom
#   type -- (str) bond type, e.g. SC1 / SC2
#   
sub echo_sidechain_bond {
	my ($i, $s, $t, $type) = @_;
	my $x = $backbone[$i] + $s;
	my $y = $backbone[$i] + $t;
	$type = "$type" . "_$aa_hash{${sequence[$i]}}";
	my $desc = sprintf("%s%d", $aa_hash{${sequence[$i]}}, $i+1);
	echo_bond($x, $y, $type, $type, $desc);
}

#
# Echo a ring bond definition.
#   i    -- (int) index of backbone atom
#   s    -- (int) ID shift of 1st atom
#   t    -- (int) ID shift of 2nd atom
#   type -- (str) bond length type, i.e. TR1 / TR2 / TR3
#
sub echo_ring_bond {
	my ($i, $s, $t, $type) = @_;
	my $x = $backbone[$i] + $s;
	my $y = $backbone[$i] + $t;
	my $length = "$type" . "_$aa_hash{${sequence[$i]}}";
	my $force = "RNG_EQL";
	my $name = sprintf("%s%d", $aa_hash{${sequence[$i]}}, $i+1);
	echo_bond($x, $y, $length, $force, $name);
}

#
# Echo a constraint definition.
#   i    -- (int) index of backbone atom
#   s    -- (int) ID shift of 1st atom
#   t    -- (int) ID shift of 2nd atom
#   type -- (str) bond length type, e.g. TR1 / TR2 / TR3
#   
sub echo_constraint {
	my ($i, $s, $t, $type) = @_;
	my $x = $backbone[$i] + $s;
	my $y = $backbone[$i] + $t;
	my $length = $parameter_hash{"BNLN_$type" . "_$aa_hash{${sequence[$i]}}"};
	my $desc = sprintf("%s%d", $aa_hash{${sequence[$i]}}, $i+1);
	printf("%d\t%d\t1\t%1.3f\t;\t%s\n", $x, $y, $length, $desc);
}

#
# Echo an angle definition.
#   i    -- (int) ID of 1st atom
#   j    -- (int) ID of 2nd atom
#   k    -- (int) ID of 3rd atom
#   type -- (str) angle/force type, e.g. SC1_GLU
#   desc -- (str) description of the angle
#   
sub echo_angle {
	my ($i, $j, $k, $type, $desc) = @_;
	printf("%d\t%d\t%d\t2\t%1.2f\t%d\t;\t%s\n", $i, $j, $k, 
		$parameter_hash{"ANGL_$type"}, $parameter_hash{"ANKB_$type"}, 
		$desc);
	if ($DEBUG) {
		print STDERR 
			"DEBUG: echo_angle() < i=$i, j=$j, k=$k, type=$type, desc=$desc\n";
	}
}

#
# Echo an improper dihedral definition.
#   i    -- (int) index of backbone atom
#   s    -- (int) ID shift of 1st atom
#   t    -- (int) ID shift of 2nd atom
#   u    -- (int) ID shift of 3rd atom
#   v    -- (int) ID shift of 4th atom
#   type -- (str) dihedral type, i.e. RG1 / RG2
#   
sub echo_improper_dihedral {
	my ($i, $s, $t, $u, $v, $type) = @_;
	my @atoms = ($backbone[$i] + $s, $backbone[$i] + $t, 
		$backbone[$i] + $u, $backbone[$i] + $v);
	$type = "$type" . "_$aa_hash{${sequence[$i]}}";
	my $desc = sprintf("%s%d", $aa_hash{${sequence[$i]}}, $i+1);
	printf("%d\t%d\t%d\t%d\t2\t%1.2f\t%d\t;\t%s\n",
		@atoms, $parameter_hash{"DIAN_$type"}, 
		$parameter_hash{"DIKB_$type"}, $desc);
}

#
# Echo a proper dihedral definition.
#   i     -- (int) index of backbone atom
#   angle -- (str) angle type, e.g. BBN_HLX
#   force -- (str) force type, e.g. BBN_HLX
#   
sub echo_proper_dihedral {
	my ($i, $angle, $force) = @_;
	my @atoms = ($backbone[$i-1], $backbone[$i], $backbone[$i+1], 
		$backbone[$i+2]);
	my $desc = sprintf("%s%d", $aa_hash{${sequence[$i]}}, $i+1);
	printf("%d\t%d\t%d\t%d\t1\t%1.2f\t%d\t1\t;\t%s\n",
		@atoms, $parameter_hash{"DIAN_$angle"}, 
		$parameter_hash{"DIKB_$force"}, $desc);
}


#---------------------------------------------------------------------------#
# MAIN LOOP
#

# default options
$ELASTIC_NETWORK = 1;
$VERBOSE = 0;
$DEBUG = 0;
$HELP = 0;
%options = ('elastic' => \$ELASTIC_NETWORK, 'help' => \$HELP, 
	'verbose' => \$VERBOSE, 'debug' => \$DEBUG);
# parse options
GetOptions(\%options, 'sequence|s=s', 'structure|ssd|2=s', 
	'topology|itp|t=s', 'custom-parameters|par|p=s', 
	'output-parameters|out|o=s', 'cysteines|cys|c=s', 'elastic|e', 
	'noelastic|dihedrals|n', 'help|h|?', 'verbose|v', 'debug|d') 
	or echo_help_hint and exit;

# help needed?
if ($HELP) {
	echo_man;
	exit;
}

if ($DEBUG) {
	print STDERR "DEBUG: options={ ";
	print STDERR map { "$_=>$options{$_}, " } keys %options;
	print STDERR " }\n";
}

# local elastic networks?
if (exists $options{noelastic}) {
	$ELASTIC_NETWORK = 0;
}

# parse arguments
$numargs = scalar @ARGV;
if ($numargs == 0) {
	if (!(exists $options{'sequence'} and exists $options{'structure'})) {
		error "Sequence and structure files have to be specified.";
	}
	$fn_seq = <$options{'sequence'}>;
	$fn_ssd = <$options{'structure'}>;
} elsif ($numargs == 2) {
	$fn_seq = <$ARGV[0]>;
	$fn_ssd = <$ARGV[1]>;
} elsif ($numargs == 1) {
	if (exists $options{'structure'}) {
		$fn_seq = <$ARGV[0]>;
		$fn_ssd = <$options{'structure'}>;
	} else {
		$fn_seq = <$options{'sequence'}>;
		$fn_ssd = <$ARGV[0]>;
	}
} else {
	error "Wrong number of arguments.";
}

# redirect STDOUT?
if (exists $options{'topology'}) {
	open STDOUT, '>', <$options{'topology'}> 
		or die "Can't redirect STDOUT: $!";
}

# be verbose about options?
if ($VERBOSE) {
	if ($ELASTIC_NETWORK) {
		print STDERR "Using local elastic networks.\n";
	} else {
		print STDERR "Using dihedrals instead of local elastic networks.\n";
	}
}


#
# READ FILES
# 

# read sequence ignoring comment lines
open SEQ, '<', $fn_seq or die "Can't open sequence file: $!";
$sequence_txt = "";
@begin_chain = ();
@end_chain = ();
$number_of_chains = 0;
while ($line = <SEQ>) {
	if ($line =~ /^>/) {
		if ($number_of_chains) {
			push @begin_chain, $b;
			$e = (length $sequence_txt) - 1;
			push @end_chain, $e;
			if ($VERBOSE) {
				print STDERR "Chain $number_of_chains has " . ($e - $b + 1) 
					. " residues (" . ($b + 1) . " - " . ($e + 1) . ").\n";
				}
		}
		$b = length $sequence_txt;
		$number_of_chains++;
	} else {
		chomp $line;
		$sequence_txt .= $line;
	}
}
# process the last chain info
push @begin_chain, $b;
$e = (length $sequence_txt) - 1;
push @end_chain, $e;
if ($VERBOSE) {
	print STDERR "Chain $number_of_chains has " . ($e - $b + 1) 
		. " residues (" . ($b + 1) . " - " . ($e + 1) . ").\n";
}
close SEQ;
# split into an array
@sequence = split(//, $sequence_txt);
$total_residues = scalar @sequence;

# read secondary structure
open SSD, '<', $fn_ssd or die "Can't open structure file: $!";
($i, $structure_txt) = <SSD>;
close SSD;
chomp $structure_txt;
# unify the syntax to use only C for coils, H for helices and E for extended
$structure_txt =~ s/~/C/g;
$structure_txt =~ s/G|I/H/g;
$structure_txt =~ s/B/E/g;
# split into an array
@structure = split(//, $structure_txt);

# check if the number of residues in the seq and ssd files match
if ($total_residues != scalar @structure) {
	error "The number of residues in the sequence file ($total_residues) " 
		. "does not match the number of residues in the structure file (" 
		. (scalar @structure) . ").";
}

# check for cysteines
if ($VERBOSE) {
	printf STDERR "Found %s cysteines.\n", scalar grep(/C/, @sequence);
}
# build cysteine bridges
if (exists $options{'cysteines'}) {
	open FILE, $options{'cysteines'} 
		or die "Can't open file " . $options{'cysteines'} . ".";
	chomp (@lines = <FILE>);
	close (FILE);
	
	if ($VERBOSE) {
		printf STDERR "I will generate %s CYS bridges.\n", scalar @lines;
	}
	@bridges = ();
	foreach $line (@lines) {
		($b, $e) = split /\s+/, $line 
			or die "Something wrong with " . $options{'cysteines'} 
			. ". Make sure there are two residues per line.\n";
		if ($sequence[$b - 1] ne "C") {
			error "Residue $b is not a cysteine. Check your .cys file.";
		}
		if ($sequence[$e - 1] ne "C") {
			error "Residue $e is not a cysteine. Check your .cys file.";
		}
		push @bridges, ($b, $e);
		if ($VERBOSE) {
			print STDERR "Making a CYS bridge between residues $b and $e.\n";
		}
	}
}

# custom parameters?
if ($options{'custom-parameters'}) {
	open FILE, $options{'custom-parameters'} 
		or die "Can't open file '" . $options{'custom-parameters'} . "'.";
	$count = 0;
	while (<FILE>) {
		chomp;
		($key, $value) = split /\s+/;
		if (exists $parameter_hash{$key}) {
			if ("$parameter_hash{$key}" ne "$value") {
				if ($VERBOSE) {
					print STDERR "Parameter $key changed from " 
						. $parameter_hash{$key} . " to $value.\n";
				}
				$parameter_hash{$key} = $value;
			}
		} else {
			error "Unknown parameter '$key' specified in " 
				. $options{'custom-parameters'} . ".";
		}
		$count++;
	}
	if ($VERBOSE) {
		print STDERR "Read $count parameters.\n";
	}
	if ($parameter_hash{ITV_BONDS} !~ /constraints|bonds/) {
		error "Parameter ITV_BONDS should be either bonds or constraints.";
	}
	if ($parameter_hash{RING_BONDS} !~ /constraints|bonds/) {
		error "Parameter RING_BONDS should be either bonds or constraints.";
	}
}

# create atoms
my $number_of_atoms=0;
@atoms = ();
for (my $i=0; $i < $total_residues; $i++) {
	$number_of_atoms += create_atoms($number_of_atoms, $i);
}
if ($VERBOSE) {
	print STDERR "Created $number_of_atoms atoms.\n";
}

# check for chains, helical stretches and switch to appropriate backbone 
# bead types
for (my $i=0; $i < $number_of_chains; $i++) {
	# parse all helices
	$number_of_helices = 0;
	@start_helix=undef;
	@end_helix=undef;
	if ($structure[$begin_chain[$i]] eq "H") {
		$number_of_helices++;
		$start_helix[$number_of_helices] = $begin_chain[$i];
		if ($structure[${begin_chain[$i]} + 1] ne "H") {
			$end_helix[$number_of_helices] = $begin_chain[$i];
		}
	}
	for (my $j=${begin_chain[$i]}+1; $j <= ${end_chain[$i]}-1; $j++) {
		if ( ($structure[$j] eq "H") && ($structure[$j-1] ne "H") ) {
			$number_of_helices++;
			$start_helix[$number_of_helices] = $j;
		}
		if ( ($structure[$j] eq "H") && ($structure[$j+1] ne "H")) {
			$end_helix[$number_of_helices] = $j;
		}
	}
	if ($structure[$end_chain[$i]] eq "H") {
		if ($structure[${end_chain[$i]}-1] eq "H") {
			$end_helix[$number_of_helices] = $end_chain[$i];
		} else {
			$number_of_helices++;
			$end_helix[$number_of_helices] = $end_chain[$i];
			$start_helix[$number_of_helices] = $end_chain[$i];;
		}
	}
	# change the bead types of peripheral helical backbones
	for (my $x=1; $x <= $number_of_helices; $x++) {
		my $s = $start_helix[$x];
		my $e = $end_helix[$x];
		my $length = $e - $s + 1;
		if ($length >= 8) {
			$atoms[${backbone[$s  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN1"};
			$atoms[${backbone[$s+1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN1"};
			$atoms[${backbone[$s+2]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN1"};
			$atoms[${backbone[$s+3]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN1"};
			$atoms[${backbone[$e-3]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN2"};
			$atoms[${backbone[$e-2]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN2"};
			$atoms[${backbone[$e-1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN2"};
			$atoms[${backbone[$e  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN2"};
		} elsif ($length == 7) {
			$atoms[${backbone[$s  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN1"};
			$atoms[${backbone[$s+1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN1"};
			$atoms[${backbone[$s+2]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN1"};
			$atoms[${backbone[$s+3]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$e-2]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN2"};
			$atoms[${backbone[$e-1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN2"};
			$atoms[${backbone[$e  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN2"};
		} elsif ($length == 6) {
			$atoms[${backbone[$s  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN1"};
			$atoms[${backbone[$s+1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN1"};
			$atoms[${backbone[$s+2]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$e-2]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$e-1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN2"};
			$atoms[${backbone[$e  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN2"};
		} elsif ($length == 5) {
			$atoms[${backbone[$s  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN1"};
			$atoms[${backbone[$s+1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$s+2]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$e-1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$e  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN2"};
		# The rest are probably redundant and insignificant and probably 
		# impossible according to dssp definitions, but for sake of 
		# completeness.....
		} elsif ($length == 4) {
			$atoms[${backbone[$s  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$s+1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$e-1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$e  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
		} elsif ($length == 3) {
			$atoms[${backbone[$s  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$s+1]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$e  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
		} elsif ($length == 2) {
			$atoms[${backbone[$s  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
			$atoms[${backbone[$e  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
		} elsif ($length == 1) {
			$atoms[${backbone[$s  ]}-1]{bead} = $parameter_hash{"BEAD_BBH_CN3"};
		}
	}
}

# change Na, Nd and Nda back to N0 in alanines and Nd and Nda to N0 in 
# prolines
for (my $i=0; $i < $number_of_atoms; $i++) {
	if ( ($atoms[$i]{resname} eq "ALA") 
		&& ($atoms[$i]{bead} =~ /Na|Nd|Nda/) ) {
		$atoms[$i]{bead} = "N0";
	}
	if ( ($atoms[$i]{resname} eq "PRO") 
		&& ($atoms[$i]{bead} =~ /Nd|Nda/) ) {
		$atoms[$i]{bead} = "N0";
	}
}

# change the termini to be charged for each chain
if ($parameter_hash{"Q_TERMINI"} eq "CHARGED") {
	for (my $j=0; $j < $number_of_chains; $j++) {
		$atoms[${backbone[$begin_chain[$j]]}-1]{bead}   = "Qd";
		$atoms[${backbone[$begin_chain[$j]]}-1]{charge} = 1.0;
		$atoms[${backbone[$end_chain[$j]  ]}-1]{bead}   = "Qa";
		$atoms[${backbone[$end_chain[$j]  ]}-1]{charge} = -1.0;
	}
} elsif ($VERBOSE) {
	print STDERR "\nCAUTION: The termini are not charged. I hope this is " 
		. "what you want.\n";
}

# update backbone bead names
foreach my $i (@backbone) {
	update_bead_name($i - 1);
}


#
# OUTPUT PARAMETERS?
#
if (exists $options{'output-parameters'}) {
	open FILE, ">", $options{'output-parameters'} 
		or die "Can't open file $options{'output-parameters'}: $!";
	foreach $key (sort keys %parameter_hash) {
		print FILE "$key\t$parameter_hash{$key}\n";
	}
	close FILE;
}


#
# OUTPUT TOPOLOGY
#
print ";;; MARTINI $martini coarse-grained topology\n";
print ";;; Generated by $my_name version $version\n\n";
print ";;; SEQ: ";
if ($VERBOSE) {
	print STDERR "\nGenerating a topology file for the following sequence: \n\n";
}
$count = 0;
for ($i=0; $i < $total_residues; $i++) {
	print $sequence[$i];
	if ($VERBOSE) { print STDERR $sequence[$i] }
	$count++;
	if ($count%10 == 0) {
		print " ";
		if ($VERBOSE) {	print STDERR " " }
	}
	if ($count%50==0) {
		print " $count\n;;; SEQ: ";
		if ($VERBOSE) { print STDERR " $count\n" }
	}
}
print "\n;;; Total Number of Amino Acid Residues: $total_residues\n";
if ($VERBOSE) {
	print STDERR "\nTotal number of amino acid residues: $total_residues\n\n";
}
print "\n\n[moleculetype]\n";
print ";molname\texclusions\n";
print "Protein\t1\n\n";

# ATOMS
print "[atoms]\n";
for my $i ( 0 .. $#atoms ) {
	print join("\t", $atoms[$i]{id}, $atoms[$i]{bead}, $atoms[$i]{resid}, 
		$atoms[$i]{resname}, $atoms[$i]{name}, $atoms[$i]{id}, 
		(sprintf "%1.3f", $atoms[$i]{charge}), ";\t$atoms[$i]{conf}\n");
}

# BONDS
print "\n[ bonds ]\n";
print ";backbone-backbone bonds\n";
for (my $j=0; $j < $number_of_chains; $j++) {
	for (my $i=$begin_chain[$j]; $i < $end_chain[$j]; $i++) {
		$curr = $ss_hash{${structure[$i]}};
		$next = $ss_hash{${structure[$i+1]}};
		$f_curr = $parameter_hash{"BNKB_BBN_$curr"};
		$f_next = $parameter_hash{"BNKB_BBN_$next"};
		# use the weakest force constant
		if ($f_curr > $f_next) {
			$force = $next;
		} else {
			$force = $curr;
		}
		echo_bond($backbone[$i], $backbone[$i+1], "BBN_$curr", 
			"BBN_$force", "$curr" . "-" . "$next");
	}
}

print ";backbone-sidechain bonds\n";
if ($parameter_hash{"ITV_BONDS"} eq "constraints") {
	for (my $i=0; $i < $total_residues; $i++) {
		if ( ($numbeads_hash{$sequence[$i]} != 1) 
			&& ($sequence[$i] !~ /I|T|V/) ) {
			echo_sidechain_bond($i, 0, 1, "SC1"); 
		}
	}
} elsif ($parameter_hash{"ITV_BONDS"} eq "bonds") {
	if ($VERBOSE) {
		print STDERR 
			"Generating bb-sc bonds instead of constraints for I/T/V.\n";
	}
	for (my $i=0; $i < $total_residues; $i++) {
		if ($numbeads_hash{$sequence[$i]} != 1) {
			echo_sidechain_bond($i, 0, 1, "SC1"); 
		}
	}
}
print ";sc1-sc2 bonds (ARG, LYS)\n";
for (my $i=0; $i < $total_residues; $i++) {
	if ($numbeads_hash{$sequence[$i]} == 3) {
		echo_sidechain_bond($i, 1, 2, "SC2"); 
	}
}

# generate bonds instead of constraints for ring structures?
if ($parameter_hash{"RING_BONDS"} eq "bonds") {
	print ";sc-sc bonds for RINGS (TRP, TYR, PHE, HIS)\n";
	if ($VERBOSE) {
		print STDERR "Generating sc-sc bonds instead of constraints for " 
			. "RINGS (TRP, TYR, PHE, HIS)\n";
	}
	for (my $i=0; $i < $total_residues; $i++) {
		if ($numbeads_hash{$sequence[$i]} > 3) {
			# HIS, PHE and TYR + TRP
			echo_ring_bond($i, 1, 2, 'TR1');
			echo_ring_bond($i, 1, 3, 'TR2');
			echo_ring_bond($i, 2, 3, 'TR3');
			if ($numbeads_hash{$sequence[$i]} > 4) {
				# TRP
				echo_ring_bond($i, 2, 4, 'TR3');
				echo_ring_bond($i, 3, 4, 'TR3');
			}
		}
	}
}

# cysteine bridges
if (exists $options{cysteines}) {
	print "; disulphide bridges (CYS-CYS)\n";
	while (@bridges) {
		my $i = shift @bridges;
		my $j = shift @bridges;
		echo_bond($backbone[$i-1] + 1, $backbone[$j-1] + 1, 
			"BRD_CYS", "BRD_CYS", sprintf("CYS%d-CYS%d", $i, $j));
	}
}

# 1-5 bonds?
if ($parameter_hash{YES_BBN_CONST15} eq "yes") {
	print ";1-5 HB constraints\n";
	for (my $j=0; $j < $number_of_chains; $j++) {
		for (my $i=$begin_chain[$j]; $i < $end_chain[$j]-3; $i++) {
			if ((substr $structure_txt, $i, 5) eq 'HHHHH') {
				my $name = $ss_hash{${structure[$i]}} . "-" 
						 . $ss_hash{${structure[$i+4]}};
				echo_bond($backbone[$i], $backbone[$i+4], "HBN_HLX", 
					"HBN_HLX", $name);
			}
		}
	}
}

# elastic network instead of dihedrals?
if ($ELASTIC_NETWORK) {
	my %short_bonds = ();
	my %long_bonds = ();
	for (my $j=0; $j < $number_of_chains; $j++) {
		for (my $i=$begin_chain[$j]+1; $i < $end_chain[$j]-1; $i++) {
			if ( substr($structure_txt, $i-1, 4) eq 'EEEE' ) {
				# atom and residue details
				($atom0, $res0) = (sprintf("%08d", $backbone[$i - 1]), 
					join('', "$aa_hash{${sequence[$i - 1]}}", $i));
				($atom1, $res1) = (sprintf("%08d", $backbone[$i]), 
					join('', "$aa_hash{${sequence[$i]}}", $i + 1));
				($atom2, $res2) = (sprintf("%08d", $backbone[$i + 1]), 
					join('', "$aa_hash{${sequence[$i + 1]}}", $i + 2));
				($atom3, $res3) = (sprintf("%08d", $backbone[$i + 2]), 
					join('', "$aa_hash{${sequence[$i + 2]}}", $i + 3));
				# store first to avoid duplicates
				$short_bonds{"$atom0 $atom2"} = "${res0}-${res2}";
				$short_bonds{"$atom1 $atom3"} = "${res1}-${res3}";
				$long_bonds{"$atom0 $atom3"} = "${res0}-${res3}";
			}
		}
	}
	# short range elastic bonds
	print ";short elastic bonds\n";
	foreach $key (sort (keys %short_bonds)) {
		my ($i, $j) = split(/ /, $key);
		echo_bond($i, $j, 'ELASTIC_SHORT', 'ELASTIC', $short_bonds{$key});
	}
	# long range elastic bonds
	print ";long elastic bonds\n";
	foreach $key (sort (keys %long_bonds)) {
		my ($i, $j) = split(/ /, $key);
		echo_bond($i, $j, 'ELASTIC_LONG', 'ELASTIC', $long_bonds{$key});
	}
}

# CONSTRAINTS
print "\n[ constraints ]\n";

if ($parameter_hash{"RING_BONDS"} eq "constraints") {
	print ";sc-sc constraints (Ring Structures)\n";
	for (my $i=0; $i < $total_residues; $i++) {
		if ($numbeads_hash{$sequence[$i]} > 3) {
			# HIS, PHE and TYR + TRP
			echo_constraint($i, 1, 2, "TR1");
			echo_constraint($i, 1, 3, "TR2");
			echo_constraint($i, 2, 3, "TR3");
			if ($numbeads_hash{$sequence[$i]} > 4) {
				# TRP
				echo_constraint($i, 2, 4, "TR3");
				echo_constraint($i, 3, 4, "TR3");
			}
		}
	}
}

if ($parameter_hash{"ITV_BONDS"} eq "constraints") {
	print ";bc-sc constraints (ITV)\n";
	for (my $i=0; $i < $total_residues; $i++) {
		if ($sequence[$i] =~ /I|T|V/) {
			echo_constraint($i, 0, 1, "SC1");
		}
	}
}

# ANGLES
print "\n[angles]\n";

# backbone angles (ss dependent)
print ";backbone-backbone-backbone angles\n";
for (my $j=0; $j < $number_of_chains; $j++) { 
	for (my $i=$begin_chain[$j]+1; $i < $end_chain[$j]; $i++) {
		my $prev = $ss_hash{${structure[$i-1]}};
		my $curr = $ss_hash{${structure[$i]}};
		my $next = $ss_hash{${structure[$i+1]}};
		my $desc = "$prev" . "-" . "$curr" . "-" . "$next";
		my $type = "";
		if ( $desc =~ /(COI|EXT|TRN|BND)-(COI|EXT|TRN|BND)-(COI|EXT|TRN|BND)/ 
				&& ($prev ne $curr || $curr ne $next) ) {
			# use coil parameters if force constants are equal
			$type = "BBN_COI";
		} elsif ($curr eq "HLX" && ${sequence[$i]} eq "P") {
			# if a proline in a helix, use special values
			$type = "BBN_HPR";
		} else {
			$type = "BBN_$curr";
		}
		echo_angle($backbone[$i-1], $backbone[$i], $backbone[$i+1], 
			$type, $desc);
	}
}

print ";backbone-backbone-sidechain angles\n";
for (my $j=0; $j < $number_of_chains; $j++) {
	for (my $i=$begin_chain[$j]; $i <= $end_chain[$j]; $i++) {
		if ($numbeads_hash{$sequence[$i]} != 1) {
			my $type = "SC1_$aa_hash{$sequence[$i]}";
			if ($i == $begin_chain[$j]) {
				my $desc = sprintf("%s-%s%d", "$ss_hash{${structure[$i+1]}}",
					"$aa_hash{${sequence[$i]}}", $i+1);
				echo_angle($backbone[$i+1], $backbone[$i], $backbone[$i]+1, 
					$type, $desc);
			} else {
				my $desc = sprintf("%s-%s%d", "$ss_hash{${structure[$i-1]}}",
					"$aa_hash{${sequence[$i]}}", $i+1);
				echo_angle($backbone[$i-1], $backbone[$i], $backbone[$i]+1, 
					$type, $desc);
			}
		}
	}
}

print ";backbone-sidechain-sidechain angles (ARG and RINGS)\n";
for (my $j=0; $j < $number_of_chains; $j++) {
	for (my $i=$begin_chain[$j]; $i <= $end_chain[$j]; $i++) {
		if ($numbeads_hash{$sequence[$i]} > 2) {
			my $type = "SC2_$aa_hash{$sequence[$i]}";
			my $desc = sprintf("%s%d", $aa_hash{${sequence[$i]}}, $i+1);
			echo_angle($backbone[$i], $backbone[$i]+1, $backbone[$i]+2, 
				$type, $desc);
			if ($numbeads_hash{$sequence[$i]} > 3 ) {
				my $type = "SC3_$aa_hash{$sequence[$i]}";
				echo_angle($backbone[$i], $backbone[$i]+1, $backbone[$i]+3, 
					$type, $desc);
			}
		}
	}
}

# DIHEDRALS
print "\n[dihedrals]\n";

# improper dihedrals (GMX Type 2, for rings)
print ";improper dihedral angles\n";
for (my $i=0; $i < $total_residues; $i++) {
	if ($numbeads_hash{$sequence[$i]} > 3) {
		echo_improper_dihedral($i, 0, 2, 3, 1, "RG1");
		if ($numbeads_hash{$sequence[$i]} > 4) {
			echo_improper_dihedral($i, 1, 2, 4, 3, "RG2");
		}
	}
}

# proper dihedrals
print ";proper dihedral angles\n";
if ($parameter_hash{YES_BBN_PROPDIH} eq "yes") {
	print ";helix backbone dihedrals\n";
	for (my $j=0; $j < $number_of_chains; $j++) {
		for (my $i=$begin_chain[$j]+1; $i < $end_chain[$j]-1; $i++) {
			if (substr($structure_txt, $i-1, 4) eq 'HHHH') {
				# if one of the middle residues is a proline
				if ( ($sequence[$i] eq "P") || ($sequence[$i+1] eq "P") ) {
					echo_proper_dihedral($i, "BBN_$ss_hash{$structure[$i]}", 
						"HLX_PRO");
				} else {
					echo_proper_dihedral($i, "BBN_$ss_hash{$structure[$i]}", 
						"BBN_$ss_hash{$structure[$i]}");
				}
			}
			# dihedrals instead of elastic networks?
			if ( !($ELASTIC_NETWORK) 
				&& (substr($structure_txt, $i-1, 4) eq 'EEEE') ) {
				echo_proper_dihedral($i, "BBN_$ss_hash{$structure[$i]}", 
					"BBN_$ss_hash{$structure[$i]}");
			}
		}
	}
}
print "\n";

if ($VERBOSE) {
	print STDERR "DONE. Topology has been generated. ";
	print STDERR "Please check it for bugs/errors...\n";
	# OBLIGATORY FUNNY LINE
	print STDERR "\n(gcq#007): \"TOPOLOGY FILE FOR MAKE BENEFIT GLORIOUS " 
		. "CG SIMULATION\"\n\n";
}

