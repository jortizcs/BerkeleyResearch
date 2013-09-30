#!/usr/bin/php
// Classify the Soda Hall SCADA tags
<?php
include ('datadump.php');
include ('populate.php');

$prefix = "/home/eecs/jortiz/scada";
$rootdir = $prefix."/UCProject_UCB_SODAHALL";

$dblink = mysql_connect('local.cs.berkeley.edu', 'acme', '410soda');
if (!$dblink )	{
	die('Could not connect to database: '.mysql_error());
}

if (!mysql_select_db("broadwin", $dblink)) {
	die('Could not connect to database: '.mysql_error());
}
echo "Connected to database successfully\n";

function addStr($list, $str) {
  foreach ($list as $l) {
    if (strcmp($l,$str) == 0) {
      //      echo "found $str\n";
      return $list;
    }
  }
  $list[] = $str;
  //  print_r($list);
  return $list;
}

$units = array(
	       'MAX_OAT'=>array('Deg F', 'Maximum outside air temp'),
	       'KWD'=>array('KW', 'kw demand'),
	       'KWH'=>array('kwh', 'consumption'),
	       'OAT'=>array('Deg F', 'outside air temp'),
	       'SKWH'=>array('cts', 'counter setpoint'),
	       'HPS'=>array('psi', 'HPS pressure'),
	       'ORH'=>array('Deg F', 'Rel. Humidity on roof'),
	       'SFM'=>array('lbs/hr', 'Steam flow meter'),
	       'CURTL'=>array('status', 'system curtailment flow'),
	       'EVENT'=>array('status', 'system event override'),
	       'HPSTM'=>array('x10lbs', 'HP Steam flow'),
	       'PRALM'=>array('pres. alarm', ''),
	       'STMON'=>array('Steam On', ''),
	       'S_S'=>array('Start/Stop', ''),
	       'OCCPY'=>array('status', 'System Occupancy Comm'),
	       'LOW_RAT'=>array('alarm', 'low rate'),
		   'LOW_RAT1'=>array('alarm', 'low rate'),
	       'LOW_RAT2'=>array('alarm', 'low rate'),
	       'LOW_RAT3'=>array('alarm', 'low rate'),
	       'LOW_RAT4'=>array('alarm', 'low rate'),
	       'SLCT_PID'=>array('Deg F', ''),
	       'SMK_ALM'=>array('alarm', 'smoke'),
		   'SMK_ALM1'=>array('alarm', 'smoke'),
	       'SMK_ALM2'=>array('alarm', 'smoke'),
	       'SMK_ALM3'=>array('alarm', 'smoke'),
	       'SMK_ALM4'=>array('alarm', 'smoke'),
	       'VAV__AVG'=>array('psi', 'avg variable air volume'),
	       'VAV__MIN'=>array('psi', 'min variable air volume'),
	       'VAV__MAX'=>array('psi', 'max variable air volume'),
	       'ART' =>array( 'Deg F', 'zone temp'),
	       'ARS'=>array('Deg F', 'zone set point'),
	       'ASO'=>array('Deg F', ''),
	       'AGN'=>array('Deg F', ''),
	       'VAV'=>array('psi', 'Var. air volume'),
	       'RVAV'=>array('psi', 'Var. air volume with reheat'),
	       'STA'=>array('Status', ''),
	       'VR'=>array('Deg F', ''),
	       'CFM'=>array('CFM', ''),
	       'CVP'=>array('psi', ''),
	       'CLV'=>array('psi', ''),
	       'DMP'=>array('psi', ''),
	       'FLT'=>array('alarm', 'fault'),
	       'H_C'=>array('economizer', ''),
	       'HVP'=>array('psi', ''),
	       'MAT'=>array('Deg F', ''),
	       'SAS'=>array('Deg F', ''),
	       'SAT'=>array('Deg F', ''),
	       'SMK'=>array('Smoke Alrm', ''),
	       'SPD'=>array('', ''),
	       'STP'=>array('w.c','unknown'),
	       'RM_SAS'=>array('Deg F','VAV driven supply air'),
	       'RAT'=>array('Deg F','Return air temp'),
	       'DP_STA'=>array('status','Differential air pressure status'),
	       'SPR'=>array('units','VFD Speed Reset'),
	       'NSB'=>array('units','Night minimum VA'),
	       'KW'=>array('kw', ''),
	       'SWS'=>array('Deg F', 'supply water temp set point'),
	       'SWT'=>array('Deg F', 'supply water temp'),
	       'TMR'=>array('min', 'delay timer'),
	       'VLV'=>array('position', 'isolation valve'),
	       'CDRWT'=>array('Deg F', 'Cond return temp'),
	       'A_M'=>array("position","auto/manual"),
	       'L_L'=>array("position","chiller lead/L"),
	       'TON'=>array("tons","primary cooling tons C"),
	       'FLOW'=>array("gpm","chilled water flow"),
	       'F_FLOW'=>array("gpm","filtered chill water flow"),
		   'D_TM'=>array("Mil/T", "Day Start Time"),
		   'N_TM'=>array("Mil/T", "Night Start Time"),
		   'BLD_EVENT'=>array("boolean","Building Event Flag"),
		   'BLD___HPS'=>array("boolean","AlarmToEmail: Email Out Test"),
		   'BLD_CURTL'=>array("boolean","Building CURTL Override"),
		   'BLD_1SKWH'=>array("Counts","Counter 1 Setpoint Counter"),
		   'BLD_STMON'=>array("boolean","Steam to Bldg on/off flag"),
		   'BLD_1_KWH'=>array("Counts","BLD kwh Meter 1"),
		   'BLD_1_KWD'=>array("KW", "KW Demand"),
		   'BLD_1_OAT'=>array("Deg F", "BLD Outside Air Temperature"),
		   'BLD_2SKWH'=>array("Counts","Counter 2 Setpoint"),
		   'BLD___SFM'=>array("Lbs/Hr", "BLD Steam Flow Meter"),
		   'BLD_HPSTM'=>array("x10Lbs", "Building HP Steam Flow To"),
		   'BLD_PRALM'=>array("boolean", "BLD Pressure Alarm"),
		   'MAX_OAT'=>array("Deg F","Max Outside Air Setpoint"),
		   'BLD_2_KWH'=>array("Counts", "BLD KWH Meter 2"),
		   'BLD_2_KWD'=>array("KW", "KW Demand"),
		   'BLD___ORH'=>array("Deg F", "RH on Roof")
);

$types = array();
$d = dir("/export/smote/local/soda-scada/UCProject_UCB_SODAHALL");
#echo "Handle: " . $d->handle . "\n";
#echo "Path: " . $d->path . "\n";
while ($entry = $d->read()) {
  $match = array();
  /*if (ereg('^SODA([0-9])R([^_]*)_+([^_].*)', $entry, $match)) {	// room stats
    $type = $match[3];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];

	if (preg_match("/^([1-9]+)/", $type, $tmatch)) {
		$type = substr($type, count($tmatch)-1, strlen($type));
		$suffix = "_".$tmatch[1];
		//echo "Wrong type, ".$tmatch[1]."\t";
	}

	$sensor_id = substr($entry, 0, strlen($entry)-strlen($type));
    echo "$sensor_id, $entry => Soda AHU: $match[1] room: $match[2] type: $type unit: $unit - $desc \n";
	
	$floor_no = floor($match[2]/100);
	$bname = "Soda";
	$type1 = "AHU";

	$sensor_id = $entry;
	$sensor_name = "Temp";
	$sensor_desc = "Temperature Sensor";
	$unit1 = "F";
	//$unit_desc = "Degrees Farenheit";
	$room_no = $match[2];
	$system_id = $match[1];

	//register($bname, $floor_no, $room_no, $type1.$system_id, $unit1, $type1, $sensor_id, $sensor_name, $sensor_desc, $unit, $unit_desc);
			
	register_zone_sensor($sensor_id, $bname, $floor_no, $room_no, $type1.$system_id, $unit, $desc, $type);

	//echo $sensor_id.": "." dumping associated data\n";
	//datadump($rootdir, $sensor_id, "Min", $dblink);
	
	
  }else if (ereg('^SODA([0-9])R([0-9]{3})(.)([A-Z].*)', $entry, $match)) { // room stats no _
    $type = $match[4];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda AHU: $match[1] room: $match[2]-$match[3] type: $type unit: $unit - $desc \n";

	$floor_no = substr($match[2],0,strlen($match[2])-2);
	$bname = "Soda";
	$type1 = "AHU";

	$sensor_id = $entry;
	$sensor_name = "Pressure";
	$sensor_desc = "Pressure Sensor";
	$unit_desc = $desc;
	$system_id = $match[1];
	$room_no = $match[2]."-".$match[3];

	//register($bname, $floor_no, $room_no, $type1.$system_id, $unit, $type1, $sensor_id, $sensor_name, $sensor_desc, $unit, $unit_desc);
	
	register_zone_sensor($sensor_id, $bname, $floor_no, $room_no, $type1.$system_id, $unit, $desc, $type);	
	
	//echo $sensor_id.": "." dumping associated data\n";
	//datadump($rootdir, $sensor_id, "Min", $dblink);

 } else if (ereg('^SODA([0-9])([0-9|A]{4})_+([A-Z].*)$', $entry, $match)) { // room stats no R
    $type = $match[3];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda AHU: $match[1] room: $match[2] type: $type unit: $unit - $desc \n";

	$system_id = $match[1];
	$room_no = $match[2];
	$floor_ = substr($match[2], 0, strlen($match[2])-1);
	$floor_no = floor($floor_/100);
	$bname = "Soda";
	$type1 = "AHU";
	$unit_desc = $desc;

	$sensor_id = $entry;
	$sensor_name = "Pressure";
	$sensor_desc = "Pressure Sensor";

	//register($bname, $floor_no, $room_no, $type1.$system_id, $unit, $type1, $sensor_id, $sensor_name, $sensor_desc, $unit, $unit_desc);
		
	register_zone_sensor($sensor_id, $bname, $floor_no, $room_no, $type1.$system_id, $unit, $desc, $type);	

	//echo $sensor_id.": "." dumping associated data\n";
	//datadump($rootdir, $sensor_id, "Min", $dblink);

  } else if (ereg('^SODA([0-9])C([^_]*)_+([^_].*)', $entry, $match)) { // corridor
    $type = $match[3];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda AHU: $match[1] corridor: $match[2] type: $type unit: $unit - $desc \n";

	$system_id = $match[1];
	$room_no = $match[2];
	$floor_no = substr($match[2], 0, 1);
	$bname = "Soda";
	$type1 = "AHU";
	$unit_desc = $desc;

	$sensor_id = $entry;
	$sensor_name = "Pressure";
	$sensor_desc = "Pressure Sensor";


	register_zone_sensor($sensor_id, $bname, $floor_no, $room_no, $type1.$system_id, $unit, $desc, $type);
  } else if (ereg('^SOD_*BLD_([0-9])_*([A-Z]*)$', $entry, $match)) { // bldg
    $type = $match[2];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda Bldg, Unit: $match[1] type: $type unit: $unit - $desc\n";

	$system_id = $match[1];
	$bname = "Soda";

	$sensor_id = $entry;

	register_zone_sensor($sensor_id, $bname, $floor_no, $room_no, $type1.$system_id, $unit, $desc, $type);
  } else if (ereg('^SOD__BLD_*([A-Z]*)$', $entry, $match)) { // bldg
    $type = $match[1];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];    
    echo "$entry => Soda Bldg, type: $type unit: $unit - $desc\n";
	////////// not included /////////

  } else if (ereg('^SOD_+(.*)$', $entry, $match)) {
    $type = $match[1];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda type: $type unit: $units - $desc\n";
	////////// not included /////////
  
  } else if (ereg('^SODA([0-9])_+(.*)$', $entry, $match)) {
    $type = $match[2];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda AHU: $match[1], type: $type unit: $unit - $desc\n";
	$system = "AHU".$match[1];
	
	register_system_sensor($entry, $system, $type);

  } else if (ereg('^SODA([0-9])S([0-9]+)_*([a-z|A-Z].*)$', $entry, $match)) { // Supply
    $type = $match[3];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda, AHU: $match[1] Supply Fan: $match[2], type: $type unit: $unit - $desc\n";
	/////////////////////////////////////////////////////
	//////////////// 95 entries	 ///////////////////////
	////////////////////////////////////////////////////
	$system = "SF".$match[2];
	$parent_sid = "AHU".$match[1];
	insert_or_update_systree($system, $parent_sid);
	register_system_sensor($entry, $system, $type);
  } else*/ if (ereg('^SODA([0-9])SF_+([A-Z].*)$', $entry, $match)) { // Supply Fan
    $type = $match[2];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda, AHU: $match[1] type: $type unit: $unit - $desc\n";
  } /*else if (ereg('^SODA([0-9])S_+([a-z|A-Z].*)$', $entry, $match)) { // Supply
    $type = $match[2];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda, Unit: $match[1] type: $type unit: $unit - $desc\n";
  } else if (ereg('^SODA([0-9])E([0-9]*)_*([a-z|A-Z].*)$', $entry, $match)) { // Exhaust
    $type = $match[3];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda, Unit: $match[1] Exhaust Fan: $match[2] type: $type unit: $unit - $desc\n";
  } else if (ereg('^SODA([0-9])ALL_*([0-9]+)([A-Z].*)$', $entry, $match)) { 
    $type = $match[3];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda, AHU: $match[1] LS$match[2] type: $type unit: $unit - $desc\n";
  } else if (ereg('^SODC([0-9])C([0-9]+)_+([A-Z].*)$', $entry, $match)) { // Cooling Chiller
    $type = $match[3];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda, Cooling $match[1] CH $match[2] system type: $type unit: $unit - $desc\n";
  } else if (ereg('^SODC([0-9])C([0-9]+)([A-Z][^_]*)+_+([A-Z].*)$', $entry, $match)) { // Cooling Chiller
    $type = $match[4];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda, Cooling $match[1] CH $match[2] $match[3] type: $type unit: $unit - $desc\n";
  } else if (ereg('^SODC([0-9])P_+([A-Z].*)$', $entry, $match)) { // Cooling Primary Chiller
    $type = $match[2];
    $unit = $units[$type][0];
    $desc  = $units[$type][1];
    echo "$entry => Soda, Cooling $match[1] Primary type: $type unit: $unit - $desc\n";


    } else if (ereg('^(SOD.*[_]+)([a-zA-Z0-9_]{1,})$', $entry, $match)) {
		$trk = preg_split("/[_]+/", $entry);
		$pieces = array();
		if (sizeof($trk) > 0) {
			for	($i=1; $i<=sizeof($trk); $i++) {
				$b= $i-1;
				$pieces[$b] = $trk[$i];
			}
			$type= implode("_", $pieces);
			$cnt = strlen($type);
			$type = substr($type, 0, $cnt-1);
		}
		//echo "nf: $entry\ttype: $type\n";
		// print_r($match);
		//     $type = $match[sizeof($match)-1];
		if($units[$type]) {
			$unit = $units[$type][0];
			$desc  = $units[$type][1];
			echo "$entry => type: $type unit: $unit - $desc\n";
		}
		else {
			echo "$entry, t_type: $type\n";
			$types = addStr($types, $type);
		}
  }*/ else {
    echo "$entry not found\n";
  }
 }
$d->close();
print_r($types);


mysql_close($dblink);
?>
