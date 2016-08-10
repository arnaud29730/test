<?php

/* php rad2deg et deg2rad */


//----------------------------------------------------------------------------------------
//ALG001 Latitude Isomérique Page 3
//PARAMETRES
//ladd latitude rad 	//pee premiere excentricte ellipsoide
//3 jeu d'essais _LatitudeIsometrique

//$latDD=0.87266462600;$e = 0.08199188998;
//print_r(_LatitudeIsometrique($latDD,$e));
//1.0055265364852

//$latDD=-0.3;$e = 0.08199188998;
//print_r(_LatitudeIsometrique($latDD,$e));
//-0.30261690063245

//$latDD=0.19998903370;$e = 0.08199188998;
//print_r(_LatitudeIsometrique($latDD,$e));
//0.20000000000894

function _LatitudeIsometrique($latDD,$pee)
{
        
        $var1 = tan( ( M_PI/4 ) + ($latDD / 2 ));
        $var3 = (1 - ($pee * ( sin($latDD) ) ) );
        $var4 = (1 + ($pee * ( sin($latDD) ) ) );
        $var5 = $var3/$var4;
        $var2 = log($var1 * pow($var5 ,( $pee / 2) ));
        return $var2;
}
//----------------------------------------------------------------------------------------
//ALG019 Projection Lambert Conique Conforme (cas tangent)
//PARAMETRES
//lodd0 longitude origine / meridien origine 	//pee premiere excentricte ellipsoide	//dga demi grand axe
// ladd0 ltitue origine
//$kfe facteur echelle
//X0,Y0 coordonnées projection point origine
//SORTIES
$loddc = 0.0;
$expprj = 0.0 ; //exposant de la projection n
$cproj = 0.0 ; //consatante de la projection c
$XS =0.0; $YS=0.0; //coordonnées pole projection

//2 jeu d'essais _LatitudeIsometrique
/*
$long0 = 0.18112808800;
$lat0 = 0.97738438100;
$k0 = 1.0;
$XOR=0.0;
$YOR=0.0;
$a = 6378388.0;
$e = 0.081991890;
*/
/*
$long0 = 0.04079234433;
$lat0 = 0.86393798000;
$k0 = 0.9998773400;
$XOR=600000.0;
$YOR=200000.0;
$a = 6378249.2;
$e = 0.0824832568;
*/

/*
_PLCC($long0,$a,$e,$lat0,$k0,$XOR,$YOR);
print_r($e);echo '<br/>';
//print_r(sin($lat0));'<br/>';
print_r($expprj);echo '<br/>';
print_r($cproj);echo '<br/>';
print_r($XS);echo '<br/>';
print_r($YS);
*/

function _PLCC( $lodd0,$dga , $pee ,$ladd0 ,$kfe , $X0 , $Y0)
{
	global $expprj,$cproj,$XS,$YS;//a passer en parametre fonction ??

	$loddc = $lodd0;
	$expprj = sin($ladd0);
//print_r($expprj);echo '<br/>';
	$cotangent = 1 / tan($ladd0);
	//$cotangent = tan(M_PI/2 - $ladd0);
	//$cotangent = tan(M_PI_2 - $ladd0);
	$expliso = exp( $expprj * _LatitudeIsometrique( $ladd0,$pee ) );
	$cproj = $kfe * _GrandeurNormale($ladd0,$dga,$pee) * $cotangent * $expliso;
	$XS = $X0;
	$YS = $Y0 + $kfe *_GrandeurNormale($ladd0,$dga,$pee) * $cotangent;

}

//----------------------------------------------------------------------------------------
//ALG021 Grandeur Normale de ellipsoide page 17/18
//PARAMETRES
//ladd latitude rad 	//pee premiere excentricte ellipsoide	//dga demi grand axe
//jeu d'essais
//$latDD=0.97738438100;$a = 6378388.0;$e = 0.081991890;
//validation
//print_r(_GrandeurNormale($latDD,$a,$e));
//obtenu = 6393174.9755385 => ok jusuq'a 4 decimale 

function _GrandeurNormale($ladd, $dga ,$pee )
{
return $dga/ sqrt( 1 - ($pee*$pee) * ( sin($ladd)*sin($ladd) ) );

}
//-----------------------------------------------------------------
//ALG003
//PARAMETRES
//pee premiere excentricte ellipsoide
//expproj n
//ctprj c constante de la projection
//loddc long origine  / meridien origine
// xsp , ysp coordone pole
//loddp , laddp coordonnes gps ws84
//SORTIE
$XL93 = 0;
$YL93 = 0;
//1 jeu d'essais 
/*
$pee= 0.0824832568;// premiere excentricte ellipsoide //e
$expproj = 0.760405966; //n
$ctprj = 11603796.9767; //c
$loddc =0.04079234433;
$xsp =600000.0;
$ysp = 5657616.674;// coordone pole
$loddp =0.14551209900;
$laddp = 0.87266462600;//coordonnes gps ws84
_Wgs84ToLambert93($pee,$expproj,$ctprj,$loddc,$xsp,$ysp,$loddp,$laddp);

print_r($XL93);echo '<br/>';//1029705.0818082
print_r($YL93);echo '<br/>';//272723.85100263
*/
function _Wgs84ToLambert93($pee,$expproj,$ctprj,$loddc,$xsp,$ysp,$loddp,$laddp)
{
	global $XL93,$YL93;
	$lisotemp = _LatitudeIsometrique( $laddp,$pee );

$var1 = $ctprj * exp( -1*$expproj*$lisotemp ); //page9
$var2 = $expproj*($loddp- $loddc);
$XL93 = $xsp +$var1 * sin( $var2) ;
$YL93 = $ysp -$var1 * cos( $var2) ;

}
//-----------------------------------------------------------------
//ALG004
//PARAMETRES
/*
$X93 = 1029705.083;
$Y93 = 272723.8490;
$pee= 0.0824832568;// premiere excentricte ellipsoide //e
$expproj = 0.760405966; //n
$ctprj = 11603796.9767; //c
$loddc =0.04079234433;//long origine  / meridien origine
$xsp=600000.0;
$ysp=5657616.674;
$epsil = pow(10,-11);//tolérance de conergence
//$epsil = pow(10,-11);//tolérance de conergence

_Lambert93ToWgs84($pee,$expproj,$ctprj,$loddc, $xsp , $ysp,$epsil,$X93,$Y93);
print_r($WGS84LON);echo 'rrrr<br/>';
print_r($WGS84LAT);echo '<bbbbbr/>';
*/

//SORTIE
$WGS84LON=0;
$WGS84LAT=0; 

function _Lambert93ToWgs84($pee,$expproj,$ctprj,$loddc, $xsp  , $ysp,$epsil,$X93,$Y93 )
{
	global $WGS84LON,$WGS84LAT;
$R = sqrt( ($X93 - $xsp )*($X93 - $xsp ) + ($Y93 - $ysp )*($Y93 - $ysp ) );
$lambda = atan2( ($X93 - $xsp ) , ($ysp - $Y93  ) );
$lontemp = $loddc + $lambda / $expproj;
//print_r($lontemp);echo '<br/>';//
//page 12
$latisotemp = -1/$expproj * log(abs($R/$ctprj));
//print_r($latisotemp);echo '<br/>';//
//ALG002 a implenete
//$latfliso = _LatitudeFLatitudeIsometrique($latisotemp,$pee,$epsil);
$WGS84LON=$lontemp;
$WGS84LAT=_LatitudeFLatitudeIsometrique($latisotemp,$pee,$epsil);;

}
//------------------------------------------------------------------------------------
//ALG002
//PARAMETRES
//$latiso
//$pee=// premiere excentricte ellipsoide //e
// $epsil //tolérance de conergence
//SORTIE


//3 jeu d'essais _LatitudeFLatitudeIsometrique
/*
print_r(_LatitudeFLatitudeIsometrique(1.00552653648,0.08199188998,pow(10,-11)));echo '<br/>';//
print_r(_LatitudeFLatitudeIsometrique(-0.30261690060,0.08199188998,pow(10,-11)));echo '<br/>';//
print_r(_LatitudeFLatitudeIsometrique(0.2,0.08199188998,pow(10,-11)));echo '<br/>';//
0.87266462599667
-0.29999999996879
0.19998903369116
ok!
*/

function _LatitudeFLatitudeIsometrique($latiso,$pee,$epsil)
{
$latfliso=0;
$la0 = 2 * atan( exp ( $latiso)) - M_PI_2;
$i = 0;
$la1_0 =2*M_PI;
$run = true;


while ($run)
	{
	$i++;
	$var1 = ( 1 + $pee*sin($la1_0) ) / ( 1 - $pee*sin($la1_0) );
	$var2 = pow($var1 , $pee / 2);
	$la1_1 = 2 * atan( $var2 * exp ($latiso) ) - M_PI_2;
	$diff = abs( $la1_1 - $la1_0);
	if ( $diff < $epsil) 
		{ $latfliso  = $la1_1;$run=false;return $latfliso;}

	$la1_0=$la1_1;
	}
}

//------------------------------
//conersvion wgs to lambert 93 / convers 3 => ok
//--------------------------
$latwgs84cnv3=deg2rad(47.7885);
$lonwgs84cnv3=deg2rad(-3.15);
$expproj = 0.7256077650; //n
$ctprj = 11754255.426; //c
$xsp = 700000.0;
$ysp = 12655612.050;
//$pee = 0.08248325676;
$pee = 0.08181919106;//wgs84
//$loddc=deg2rad(2+(20/60)+(14.025/3600));
$loddc=deg2rad(3);
$epsil = pow(10,-11);//tolérance de conergence


_Wgs84ToLambert93($pee,$expproj,$ctprj,$loddc,$xsp,$ysp,$lonwgs84cnv3,$latwgs84cnv3);

print_r($XL93);echo '<br/>';//1029705.0818082
print_r($YL93);echo '<br/>';//272723.85100263
//---------------------------------------
//reciporque



_Lambert93ToWgs84($pee,$expproj,$ctprj,$loddc, $xsp , $ysp,$epsil,$XL93,$YL93);
print_r(rad2deg($WGS84LON));echo '<br/>';
print_r(rad2deg($WGS84LAT));echo '<br/>';
/*
_Lambert93ToWgs84($pee,$expproj,$ctprj,$loddc, $xsp , $ysp,$epsil,$X93,$Y93);
print_r($WGS84LON);echo 'rrrr<br/>';
print_r($WGS84LAT);echo '<bbbbbr/>';
*/
?>
