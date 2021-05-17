<b>Dictionary for Supplemental Tables S4-S8</b>
<em>(Run and table column descriptions and dictionary)</em>

<table>						
<thead><tr><th>	FileName	</th><th>	Description	</th></tr></thead>		
<tbody>						
<tr><td> 	Supplemental Table S4.csv	</td><td>	Evaluated results of Monte Carlo runs on each scenario in each year, statewide.	</td></tr>		
<tr><td> 	Supplemental Table S5.csv	</td><td>	Evaluated results of Monte Carlo runs on each pathway activity per scenario in each year, statewide.	</td></tr>		
<tr><td> 	Supplemental Table S6.csv	</td><td>	Evaluated results of Monte Carlo runs on each scenario in each year, per county.	</td></tr>		
<tr><td> 	Supplemental Table S7.csv	</td><td>	Evaluated results of Monte Carlo runs on each pathway subactivity per scenario in each year, per county.	</td></tr>		
<tr><td> 	Supplemental Table S8.csv	</td><td>	Evaluated results of Monte Carlo runs on each pathway activity per scenario in each year, per county.	</td></tr>		
</tbody>						
</table>						
<table>						
<thead><tr><th>	Columns	</th><th>	Description	</th><th>	Values or Units	</th></tr></thead>
<tbody>						
<tr><td> 	Scenario	</td><td>	Numeric identifier for each scenario.	</td><td>	2=Ambitious, 3=Moderate, 1=Limited	</td></tr>
<tr><td> 	ScenName	</td><td>	Name of each scenario	</td><td>	Amb=Ambitious, Mod=Moderate, Lim=Limited	</td></tr>
<tr><td> 	County	</td><td>	Name of each county (39 counties total)	</td><td>	Adams, Asotin, Benton, Chelan, Clallam, Clark, Columbia, Cowlitz, Douglas, Ferry, Franklin, Garfield, Grant, Grays Harbor, Island, Jefferson, King, Kitsap, Kittitas, Klickitat, Lewis, Lincoln, Mason, Okanogan, Pacific, Pend Oreille, Pierce, San Juan, Skagit, Skamania, Snohomish, Spokane, Stevens, Thurston, Wahkiakum, Walla Walla, Whatcom, Whitman, Yakima	</td></tr>
<tr><td> 	Activity	</td><td>	NCS pathway, where some are a combination of subactivities (defined below)	</td><td>	ACF=avoided conversion of forest, ag=agricultural activities, grass=restoration and avoided conversion of grasslands, ifm=improved forest management (aka extended timber harvest rotations), reforest=replanting after wildfires, rref=riparian reforestation, sage=restoration and avoided conversion of sagebrush-steppe, tide=tidal wetland restoration	</td></tr>
<tr><td> 	Subactivity	</td><td>	Specific NCS activity that is equivalent to or combines with other subactivities to be a NCS pathway	</td><td>	(see Lines 19-37 below)	</td></tr>
<tr><td> 	Year	</td><td>	Consecutive years that the pathway is implemented per each scenario	</td><td>	Range: 1-31 (e.g. representing hypothetical years 2020-2050)	</td></tr>
<tr><td> 	Median	</td><td>	Median of Monte Carlo results in MMT CO2e per year.	</td><td>	Units: MMT CO2e yr-1	</td></tr>
<tr><td> 	CI_High	</td><td>	Confidence Interval at 95th percentile value of Monte Carlo results.	</td><td>	Units: MMT CO2e yr-1	</td></tr>
<tr><td> 	CI_Low	</td><td>	Confidence Interval at 5th percentile value of Monte Carlo results.	</td><td>	Units: MMT CO2e yr-1	</td></tr>
</tbody>						
</table>						
<table>						
<thead><tr><th>	Subactivity (Line 13 above)	</th><th>	Activity (Line 12 above)	</th><th>	Description	</th></tr></thead>
<tbody>						
<tr><td> 	1f2rur	</td><td>	ACF	</td><td>	avoided conversion of forests to rural development	</td></tr>
<tr><td> 	2f2urb	</td><td>	ACF	</td><td>	avoided conversion of forests to urban development	</td></tr>
<tr><td> 	1cc	</td><td>	ag	</td><td>	cover crop	</td></tr>
<tr><td> 	2nt	</td><td>	ag	</td><td>	no-till	</td></tr>
<tr><td> 	3nm	</td><td>	ag	</td><td>	nutrient management	</td></tr>
<tr><td> 	1grass	</td><td>	grass	</td><td>	restoration and avoided conversion of grassland	</td></tr>
<tr><td> 	1feder	</td><td>	ifm	</td><td>	extended timber harvest on federal lands	</td></tr>
<tr><td> 	2other	</td><td>	ifm	</td><td>	extended timber harvest on other (e.g. nature preserves, local government) lands	</td></tr>
<tr><td> 	3priva	</td><td>	ifm	</td><td>	extended timber harvest on private lands (industrial and non-industrial)	</td></tr>
<tr><td> 	4state	</td><td>	ifm	</td><td>	extended timber harvest on state lands	</td></tr>
<tr><td> 	5prseq	</td><td>	ifm	</td><td>	extended timber harvest on private lands, with additional sequestration by avoiding clearcutting	</td></tr>
<tr><td> 	1intensL	</td><td>	reforest	</td><td>	replanting forest after low intensity wildfire	</td></tr>
<tr><td> 	2intensM	</td><td>	reforest	</td><td>	replanting forest after medium intensity wildfire	</td></tr>
<tr><td> 	3intensH	</td><td>	reforest	</td><td>	replanting forest after high intensity wildfire	</td></tr>
<tr><td> 	1rref	</td><td>	rref	</td><td>	riparian reforestation	</td></tr>
<tr><td> 	1s2ag	</td><td>	sage	</td><td>	restoration of sagebrush-steppe from agricultural uses	</td></tr>
<tr><td> 	2ag2s	</td><td>	sage	</td><td>	avoided conversation of sagebrush-steppe to agricultural uses	</td></tr>
<tr><td> 	1tide	</td><td>	tide	</td><td>	tidal wetland restoration	</td></tr>
</tbody>						
</table>						
