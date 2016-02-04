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
my $str = '2 chain types (1,2) with solvent (3)';
my $N=100;    #length of the chains
my $M=31;     #number of chains
my $R=6250;     #number of solvent
my $C=1.76;
my $b=0.97;    #bond length
my $len=25;    #cell size
my (%chains,%sys);
my $type1=1;
my $type2=1;

################## make walk - chains generation ###################
sub MakeWalk(%){
#making system hash
my %chains= ('header' => $str,
            'N' => $N,
            'M' => $M,
            'R' => $R,
            'Cinf' => $C,
            'b' => $b,
            'len' => $len);

my ($i,$k,$Q,$fi,$l,@randpos);
my @nullpos=(0,0,0);

#calculating of critical parameters
my $O_max=2*acos(sqrt(($chains{'Cinf'}-1)/($chains{'Cinf'}+1)));
my $l_max=$chains{'b'}*sqrt((1-cos(pi-$O_max))**2+sin($O_max)**2);

#make first atom
for ($i=1; $i<=$chains{'M'}; $i++) {
    @{$chains{'mols'}[$i][1]} = @nullpos;
    }
    
#make first atoms of solvent
for ($i=($chains{'M'}+1); $i<=($chains{'M'}+$chains{'R'}); $i++) {
    @{$chains{'mols'}[$i][1]} = @nullpos;
    }

#make second atom
for ($i=1; $i<=$chains{'M'}; $i++) {
    $fi = 2*pi*rand();
    $Q = pi*rand();
    @randpos = ($chains{'b'}*sin($Q)*cos($fi),$chains{'b'}*sin($Q)*sin($fi),$chains{'b'}*cos($Q));
    @{$chains{'mols'}[$i][2]} = @randpos;
    }
#make atoms from 3 to N
for ($k=3; $k<=$chains{'N'}; $k++){
    for ($i=1; $i<=$chains{'M'}; $i++) {
        my $fuse=0;
           do{
              $fuse++;

              $fi = 2*pi*rand();
              $Q = pi*rand();
                                  $randpos[0] = $chains{'mols'}[$i][$k-1][0] + $chains{'b'}*sin($Q)*cos($fi);
                                  $randpos[1] = $chains{'mols'}[$i][$k-1][1] + $chains{'b'}*sin($Q)*sin($fi);
                                  $randpos[2] = $chains{'mols'}[$i][$k-1][2] + $chains{'b'}*cos($Q);
                                 $l=sqrt(($chains{'mols'}[$i][$k-2][0]-$randpos[0])**2+($chains{'mols'}[$i][$k-2][1]-$randpos[1])**2+($chains{'mols'}[$i][$k-2][2]-$randpos[2])**2);

                                 if ($fuse>1000) { print "fuse is over N= $k M= $i"; last;}
                                 if ($l>1.94) {print "ERROR";}
              } until ($l>$l_max) ;
           @{$chains{'mols'}[$i][$k]} = @randpos;
           }
       }
return %chains;
}
############# end make walk #################
#print "Scitical params = ", rad2deg($O_max)," ", $l_max,"\n";
#################  output in console  ###################################################################

# %chains transmited as link don't work with solvent
sub out_in_console{
  my %inhash = %{shift()};
  my ($i,$k);
     for ($i=1; $i<=$inhash{'M'}; $i++){
         my $tempR2 = ($inhash{'mols'}[$i][1][0]-$inhash{'mols'}[$i][$N][0])**2+($inhash{'mols'}[$i][1][1]-$inhash{'mols'}[$i][$N][1])**2+($inhash{'mols'}[$i][1][2]-$inhash{'mols'}[$i][$N][2])**2;
         print "Molecular $i R**2 = $tempR2 \n";
               for ($k=1;$k<=$inhash{'N'};$k++){
               my @temp=@{$inhash{mols}[$i][$k]};
               print "@temp \n";
               }

}
}
#################### end output in console ###############################################################



############# output to file ######################### don't work with solvent
sub output_to_file($$){

  my %inhash = %{shift()};
  my $filename = shift();

open (my $ChainDump, '>', $filename) or die "Can't open dump file $filename \n";
print $ChainDump "This is RANDOMWalk dump file \n";


 for (my $i=1; $i<=$inhash{'M'}; $i++){
         my $tempR2 = ($inhash{'mols'}[$i][1][0]-$inhash{'mols'}[$i][$N][0])**2+($inhash{'mols'}[$i][1][1]-$inhash{'mols'}[$i][$N][1])**2+($inhash{'mols'}[$i][1][2]-$inhash{'mols'}[$i][$N][2])**2;
         print $ChainDump "Molecular $i R**2 = $tempR2 \n";
               for (my $k=1;$k<=$inhash{'N'};$k++){
               my @temp=@{$inhash{mols}[$i][$k]};
               print $ChainDump "@temp \n";
               }
 }

close $ChainDump;
}
############# end output to file #########################

############## output to pdb ######################
sub output_to_pdb($$){

  my %inhash = %{shift()};
  my $filename = shift();

open (my $ChainDump, '>', $filename) or die "Can't open dump file $filename \n";
print $ChainDump "TITLE This is RANDOMWalk dump file \n";
print $ChainDump "REMARK generated dy RANDOMWalk\n";
print $ChainDump "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n";
print $ChainDump "MODEL     1 \n";


 for (my $i=1; $i<=$inhash{'M'}; $i++){
         my $tempR2 = ($inhash{'mols'}[$i][1][0]-$inhash{'mols'}[$i][$N][0])**2+($inhash{'mols'}[$i][1][1]-$inhash{'mols'}[$i][$N][1])**2+($inhash{'mols'}[$i][1][2]-$inhash{'mols'}[$i][$N][2])**2;
         #print $ChainDump "Molecular $i R**2 = $tempR2 \n";
               for (my $k=1;$k<=$inhash{'N'};$k++){
               my @temp=@{$inhash{mols}[$i][$k]};
               my $name=1;
               my $res='MOL';
               my $resNum=1;
               printf $ChainDump "ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f". "  1.00  0.00\n", $k, $name, $res, $resNum,
                $inhash{mols}[$i][$k][0], $inhash{mols}[$i][$k][1], $inhash{mols}[$i][$k][2];
               }
 }
for (my $k=1;$k<=$inhash{'N'};$k++){

  for ($k<$inhash{'N'} and $k>1) { print $ChainDump "CONECT $k ",($k+1)," ",($k-1)," \n";}
  #else if($k=1) print $ChainDump "CONECT $k ()";

}
print $ChainDump "TER \nENDMDL \n";
close $ChainDump;
}
#####################################


##########
sub location_to_cell($){
my %chains = %{shift()};
# randomising starting position for polymer 
for (my $i=1; $i<=$chains{'M'}; $i++){
  my ($x0,$y0,$z0);
  $x0=$chains{'len'}*rand();
  $y0=$chains{'len'}*rand();
  $z0=$chains{'len'}*rand();
               for (my $k=1;$k<=$chains{'N'};$k++){
                   
                    $chains{'mols'}[$i][$k][0]+=$x0;
                    $chains{'mols'}[$i][$k][1]+=$y0;
                    $chains{'mols'}[$i][$k][2]+=$z0;
               }
}

# randomising starting position for solvent
for (my $i=($chains{'M'}+1); $i<=($chains{'M'}+$chains{'R'}); $i++){
  my ($x0,$y0,$z0);
                   $x0=$chains{'len'}*rand();
                   $y0=$chains{'len'}*rand();
                   $z0=$chains{'len'}*rand();
                   $chains{'mols'}[$i][1][0]+=$x0;
                   $chains{'mols'}[$i][1][1]+=$y0;
                   $chains{'mols'}[$i][1][2]+=$z0;
}

############ adding periodic boundary condition x,y,z=len
for (my $i=1; $i<=$chains{'M'}; $i++){
               for (my $k=1;$k<=$chains{'N'};$k++){

                   if($chains{'mols'}[$i][$k][0] > $chains{'len'}) {
                   for (my $j=$k;$j<=$chains{'N'};$j++) {$chains{'mols'}[$i][$j][0] -= $len;}
                   }
                   if ($chains{'mols'}[$i][$k][0] < 0){
                   for (my $j=$k;$j<=$chains{'N'};$j++) {$chains{'mols'}[$i][$j][0] += $len;}
                   }
                   if ($chains{'mols'}[$i][$k][1] > $chains{'len'}) {
                   for (my $j=$k;$j<=$chains{'N'};$j++) {$chains{'mols'}[$i][$j][1] -= $len;}
                   }
                   if ($chains{'mols'}[$i][$k][1] < 0){
                   for (my $j=$k;$j<=$chains{'N'};$j++) {$chains{'mols'}[$i][$j][1] += $len;}
                   }
                   if ($chains{'mols'}[$i][$k][2] > $chains{'len'}) {
                   for (my $j=$k;$j<=$chains{'N'};$j++) {$chains{'mols'}[$i][$j][2] -= $len;}
                   }
                   if ($chains{'mols'}[$i][$k][2] < 0){
                   for (my $j=$k;$j<=$chains{'N'};$j++) {$chains{'mols'}[$i][$j][2] += $len;}
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
    printf FILE "%d atoms\n", $chains{'M'}*$chains{'N'}+$chains{R};
    printf FILE "%d bonds\n", $chains{'M'}*($chains{'N'}-1);

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

    # Atoms
        printf FILE "Atoms\n\n";
        for (my $i=1; $i<=$chains{'M'}; $i++){
               for (my $k=1; $k<=$chains{'N'}; $k++){
                my $AtomNumber=$k+$chains{'N'}*($i-1);
                my $MolNumber = $i;
                #printf File $AtomNumber,"\n";
                if ($i % ($type1 + $type2) < $type1) {$type = 1;}
                else {$type = 2;}

                @temp = @{$chains{'mols'}[$i][$k]};
                 printf FILE " $AtomNumber $MolNumber $type $chains{'mols'}[$i][$k][0] $chains{'mols'}[$i][$k][1] $chains{'mols'}[$i][$k][2]\n";
        }
        }
        #writing solvent
        for (my $i=($chains{'M'}+1); $i<=($chains{'M'}+$chains{'R'}); $i++){
                my $k=1;
                my $AtomNumber=($i-$chains{'M'})+$chains{'N'}*$chains{'M'};
                my $MolNumber = $i;
                $type = 3;
                @temp = @{$chains{'mols'}[$i][$k]};
                 printf FILE " $AtomNumber $MolNumber $type $chains{'mols'}[$i][$k][0] $chains{'mols'}[$i][$k][1] $chains{'mols'}[$i][$k][2]\n";

        }
        printf FILE "\n";

    # Bonds
        printf FILE "Bonds\n\n";


        for (my $i=1; $i<=$chains{'M'}; $i++){
               for (my $k=1;$k<$chains{'N'};$k++){
                my $AtomNumber=$k+$chains{'N'}*($i-1);
                my $MolNumber = $i;

                if ($i % ($type1 + $type2) < $type1) {$type = 1;}
                else {$type = 2;}

                 my $s = $k+($chains{'N'}-1)*($i-1);
                 my $f1 = $k + ($i-1)*($chains{'N'});
                 my $f2 = $k + ($i-1)*($chains{'N'})+1 ;
                 print FILE " $s $type $f1 $f2\n";
        }
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
%chains=MakeWalk();
#writeLammpsTypeBond(\%chains, 'chainsbefore.lmps');
#convert_lammps_to_pdb('chainsbefore.lmps', 'chainsbefore.pdb');
#convert_lammps_to_psf('chainsbefore.lmps', 'chainsbefore.psf');
#out_in_console(\%chains);
#output_to_file (\%chains , 'dump.txt');
#output_to_pdb (\%chains , 'dump.pdb');
location_to_cell(\%chains);
#output_to_file (\%chains , 'dumpPBCX.txt');
writeLammpsTypeBond(\%chains, 'chains_small.lmps');
convert_lammps_to_pdb('chains_small.lmps', 'chains_small.pdb');
convert_lammps_to_psf('chains_small.lmps', 'chains_small.psf');
###############  end main ###########################

