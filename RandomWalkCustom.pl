#2 chain types (1,2) with solvent (3)
#
use strict;
use lib 'C:\Strawberry\perl\lib';
use ModPolymatic();

use Getopt::Long();
use Math::Trig;
# Import constants pi2, pip2, pip4 (2*pi, pi/2, pi/4).
use Math::Trig ':pi';

######################## system params ###############
my $str = '#2 chain types (1,2) with solvent (3) dense system';
#my $N=100;    #length of the chains
#my $M=240;     #number of chains
#my $R=2700;     #number of solvent
my $C=1.76;
my $b=0.97;    #bond length
my $len=5;    #cell size
my (%chains,%sys);
my $type1=1;
my $type2=1;

################## make walk - chains generation ###################
sub MakeWalk(%){
  #MakeWalk(M,N,b,C,Type)
my $M = shift();
my $N = shift();
my $b = shift();
my $C = shift();
my $Type = shift();
#making system hash
my %chains= ('header' => $str,
            'M' => 0,
            'R' => 0,
            'len' => 0);

my ($i,$k,$Q,$fi,$l,@randpos);
# (Type,X,Y,Z)
my @nullpos=($Type,0,0,0);

#calculating of critical parameters
my $O_max=2*acos(sqrt(($C-1)/($C+1)));
my $l_max=$b*sqrt((1-cos(pi-$O_max))**2+sin($O_max)**2);

#make first atom
for ($i=($chains{'M'}+1); $i<=($chains{'M'}+$M); $i++) {
    @{$chains{'mols'}[$i][1]} = @nullpos;
    }
if ($N > 1) {
#make second atom
for ($i=($chains{'M'}+1); $i<=($chains{'M'}+$M); $i++) {
    $fi = 2*pi*rand();
    $Q = pi*rand();
    @randpos = ($Type, $b*sin($Q)*cos($fi),$b*sin($Q)*sin($fi),$b*cos($Q));
    @{$chains{'mols'}[$i][2]} = @randpos;
    }
    if($N > 2){
#make atoms from 3 to N
for ($k=3; $k<=$N; $k++){
    for ($i=($chains{'M'}+1); $i<=($chains{'M'}+$M); $i++) {
        my $fuse=0;
           do{
              $fuse++;
              $fi = 2*pi*rand();
              $Q = pi*rand();
              $randpos[0] = $Type;
              $randpos[1] = $chains{'mols'}[$i][$k-1][1] + $b*sin($Q)*cos($fi);
              $randpos[2] = $chains{'mols'}[$i][$k-1][2] + $b*sin($Q)*sin($fi);
              $randpos[3] = $chains{'mols'}[$i][$k-1][3] + $b*cos($Q);
              $l=sqrt(($chains{'mols'}[$i][$k-2][1]-$randpos[1])**2+($chains{'mols'}[$i][$k-2][2]-$randpos[2])**2+($chains{'mols'}[$i][$k-2][3]-$randpos[3])**2);
              if ($fuse>1000) { print "fuse is over N= $k M= $i"; print $l; last;}
              if ($l>1.94) {print "ERROR too big lenght of bond";}
              } until ($l>$l_max) ;
           @{$chains{'mols'}[$i][$k]} = @randpos;
           }
       }
}
}
#update chains info
$chains{'M'} += $M;
$chains{'R'} += $N*$M;
return %chains;
}
############# end make walk #################
#AddWalk(%chains,M,N,b,C,Type)
sub AddWalk(%){
my ($chains, $M, $N, $b, $C, $Type) = @_;
my ($i,$k,$Q,$fi,$l,@randpos);
# (Type,X,Y,Z)
my @nullpos=($Type,0,0,0);

#calculating of critical parameters
my $O_max=2*acos(sqrt(($C-1)/($C+1)));
my $l_max=$b*sqrt((1-cos(pi-$O_max))**2+sin($O_max)**2);

#make first atom
for ($i=($chains{'M'}+1); $i<=($chains{'M'}+$M); $i++) {
    @{$chains{'mols'}[$i][1]} = @nullpos;
    }
if ($N > 1) {
#make second atom
for ($i=($chains{'M'}+1); $i<=($chains{'M'}+$M); $i++) {
    $fi = 2*pi*rand();
    $Q = pi*rand();
    @randpos = ($Type, $b*sin($Q)*cos($fi),$b*sin($Q)*sin($fi),$b*cos($Q));
    @{$chains{'mols'}[$i][2]} = @randpos;
    }
    if($N > 2){
#make atoms from 3 to N
for ($k=3; $k<=$N; $k++){
    for ($i=($chains{'M'}+1); $i<=($chains{'M'}+$M); $i++) {
        my $fuse=0;
           do{
              $fuse++;
              $fi = 2*pi*rand();
              $Q = pi*rand();
              $randpos[0] = $Type;
              $randpos[1] = $chains{'mols'}[$i][$k-1][1] + $b*sin($Q)*cos($fi);
              $randpos[2] = $chains{'mols'}[$i][$k-1][2] + $b*sin($Q)*sin($fi);
              $randpos[3] = $chains{'mols'}[$i][$k-1][3] + $b*cos($Q);
              $l=sqrt(($chains{'mols'}[$i][$k-2][1]-$randpos[1])**2+($chains{'mols'}[$i][$k-2][2]-$randpos[2])**2+($chains{'mols'}[$i][$k-2][3]-$randpos[3])**2);
              if ($fuse>1000) { print "fuse is over N= $k M= $i"; print $l; last;}
              if ($l>1.94) {print "ERROR too big lenght of bond";}
              } until ($l>$l_max) ;
           @{$chains{'mols'}[$i][$k]} = @randpos;
           }
       }
}
}
#update chains info
$chains{'M'} += $M;
$chains{'R'} += $N*$M;
}

##########
sub location_to_cell{
my ($chains, $len) = @_;
$chains{'len'} = $len;
# randomising starting position for polymer
for (my $i=1; $i<=$chains{'M'}; $i++){
  my ($x0,$y0,$z0);
  $x0=$chains{'len'}*rand();
  $y0=$chains{'len'}*rand();
  $z0=$chains{'len'}*rand();
               for (my $k=1;$k<=$#{$chains{'mols'}[$i]};$k++){
                   
                    $chains{'mols'}[$i][$k][1]+=$x0;
                    $chains{'mols'}[$i][$k][2]+=$y0;
                    $chains{'mols'}[$i][$k][3]+=$z0;
               }
}

############ adding periodic boundary condition x,y,z=len
for (my $i=1; $i<=$chains{'M'}; $i++){
               for (my $k=1;$k<=$#{$chains{'mols'}[$i]};$k++){

                   if($chains{'mols'}[$i][$k][1] > $chains{'len'}) {
                   for (my $j=$k;$j<=$#{$chains{'mols'}[$i]};$j++) {$chains{'mols'}[$i][$j][1] -= $len;}
                   }
                   if ($chains{'mols'}[$i][$k][1] < 0){
                   for (my $j=$k;$j<=$#{$chains{'mols'}[$i]};$j++) {$chains{'mols'}[$i][$j][1] += $len;}
                   }
                   if ($chains{'mols'}[$i][$k][2] > $chains{'len'}) {
                   for (my $j=$k;$j<=$#{$chains{'mols'}[$i]};$j++) {$chains{'mols'}[$i][$j][2] -= $len;}
                   }
                   if ($chains{'mols'}[$i][$k][2] < 0){
                   for (my $j=$k;$j<=$#{$chains{'mols'}[$i]};$j++) {$chains{'mols'}[$i][$j][2] += $len;}
                   }
                   if ($chains{'mols'}[$i][$k][3] > $chains{'len'}) {
                   for (my $j=$k;$j<=$#{$chains{'mols'}[$i]};$j++) {$chains{'mols'}[$i][$j][3] -= $len;}
                   }
                   if ($chains{'mols'}[$i][$k][3] < 0){
                   for (my $j=$k;$j<=$#{$chains{'mols'}[$i]};$j++) {$chains{'mols'}[$i][$j][3] += $len;}
                   }
}
}
}
#############  adding periodic boundary condition end





# writeLammps( \%sys, $file )
# Write LAMMPS data file for given molecular system
sub writeLammpsTypeBond
{
    # Variables
    my %chains = %{shift()};
    my $file = shift();

    my ($num, $mol, $type, $q, @temp);


    # Open file
    open FILE, "> $file" or die "Error opening file '$file': $!";

    # Header
    printf FILE "%s\n\n", $chains{'header'};

    # Counts
    printf FILE "%d atoms\n", $chains{'R'};
    printf FILE "%d bonds\n", $chains{'R'}-$chains{'M'};

    printf FILE "%d atom types\n", 3;
    printf FILE "%d bond types\n", 2;
    printf FILE "\n";

    # Box dimensions
    printf FILE "%f %f xlo xhi\n",0,$chains{'len'};
    printf FILE "%f %f ylo yhi\n",0,$chains{'len'};
    printf FILE "%f %f zlo zhi\n\n",0,$chains{'len'};

    # Masses
        printf FILE "Masses\n\n";
        printf FILE " 1 1\n 2 1\n 3 1\n ";
        printf FILE "\n";

    # Atoms       $#{$chains{'mols'}[1]} atoms in molecula
        printf FILE "Atoms\n\n";
        my  $AtomNumber = 0;
        for (my $i=1; $i<=$chains{'M'}; $i++){
               for (my $k=1; $k<=$#{$chains{'mols'}[$i]}; $k++){
                $AtomNumber++;
                my $MolNumber = $i;
                #printf File $AtomNumber,"\n";
                @temp = @{$chains{'mols'}[$i][$k]};
                 printf FILE " $AtomNumber $MolNumber $chains{'mols'}[$i][$k][0] $chains{'mols'}[$i][$k][1] $chains{'mols'}[$i][$k][2] $chains{'mols'}[$i][$k][3]\n";
        }
        }
        printf FILE "\n";

    # Bonds
        printf FILE "Bonds\n\n";
        my  $AtomNumber = 0;
        my  $BondNumber = 0;
        for (my $i=1; $i<=$chains{'M'}; $i++){
               for (my $k=1;$k<$#{$chains{'mols'}[$i]};$k++){
                $BondNumber++;
                $AtomNumber++;
                my $MolNumber = $i;

                 my $s = $BondNumber;
                 my $f1 = $AtomNumber;
                 my $f2 = $AtomNumber+1 ;
                 printf FILE " $s $chains{'mols'}[$i][$k][0] $f1 $f2\n";
        }
        #atom without bond
        $AtomNumber++;
        }
        printf FILE "\n";



    # Close file
    close FILE;
}

sub convert_lammps_to_pdb{
my %sys = Polymatic::readLammpsTypeBond($_[0]);
Polymatic::writePdb($_[1], \%sys);

}
sub convert_lammps_to_psf{
my %sys = Polymatic::readLammpsTypeBond($_[0]);
Polymatic::writePsf($_[1], \%sys);

}
###############   main  #############################
%chains=MakeWalk(3,20,$b,$C,1);
AddWalk(\%chains,3,2,$b,$C,2);
writeLammpsTypeBond(\%chains, 'chainsbefore.lmps');
convert_lammps_to_pdb('chainsbefore.lmps', 'chainsbefore.pdb');
#convert_lammps_to_psf('chainsbefore.lmps', 'chainsbefore.psf');
location_to_cell(\%chains,$len);
writeLammpsTypeBond(\%chains, 'chains_small.lmps');
convert_lammps_to_pdb('chains_small.lmps', 'chains_small.pdb');
#convert_lammps_to_psf('chains_small.lmps', 'chains_small.psf');
###############  end main ###########################

