<map version="0.8.1">
<!-- To view this file, download free mind mapping software FreeMind from http://freemind.sourceforge.net -->
<node CREATED="1226224632383" ID="Freemind_Link_1526716914" MODIFIED="1226225213899" TEXT="metmask">
<node CREATED="1226225699828" ID="Freemind_Link_1199700841" MODIFIED="1226225837846" POSITION="right" TEXT="Current method is Flawed" VSHIFT="4">
<node CREATED="1226225845918" ID="Freemind_Link_1658980022" MODIFIED="1226225886950" TEXT="Search MMDC with whatever we have and then manually curate"/>
<node CREATED="1226225887866" ID="Freemind_Link_630779839" MODIFIED="1226225896458" TEXT="slow and error-prone"/>
<node CREATED="1226225897105" ID="Freemind_Link_1202886987" MODIFIED="1226225918981" TEXT="cant edit MMDC so lots of &quot;manual fixes&quot;"/>
<node CREATED="1226225922365" ID="Freemind_Link_1759763777" MODIFIED="1226225932445" TEXT="no influence on the future of mmdc"/>
<node CREATED="1226225932949" ID="Freemind_Link_1210352539" MODIFIED="1226225940876" TEXT="cant integrate with whatever databases we have here"/>
</node>
<node CREATED="1226224680494" ID="_" MODIFIED="1226224682855" POSITION="right" TEXT="Purpose">
<node CREATED="1226224686918" ID="Freemind_Link_589259681" MODIFIED="1226224814221" TEXT="Quick retrieval of bio met from weak ids">
<node COLOR="#338800" CREATED="1226224820333" ID="Freemind_Link_1153754577" MODIFIED="1226224904485" TEXT="a weak identifier is something that can not be easily mapped to a bio met, ie string or multiple ambigous identifiers"/>
<node CREATED="1231827534425" ID="Freemind_Link_995681494" MODIFIED="1231827552398" TEXT="usefulness of weak identifers is questionable"/>
</node>
<node CREATED="1226225971046" ID="Freemind_Link_1520310429" MODIFIED="1226226003662" TEXT="not a search tool for metabolites, but for metabolite IDs"/>
<node CREATED="1226224940570" ID="Freemind_Link_1238062658" MODIFIED="1226224961009" TEXT="tidy source for all bio mets">
<node CREATED="1226225955283" ID="Freemind_Link_965500898" MODIFIED="1226225968586" TEXT="Integrate with all our current dbs and stanard lists"/>
</node>
</node>
<node CREATED="1226224998151" ID="Freemind_Link_1553056211" MODIFIED="1226225000653" POSITION="right" TEXT="components">
<node CREATED="1226225004248" ID="Freemind_Link_1200140184" MODIFIED="1226225008081" TEXT="Analytes">
<node CREATED="1226225028825" ID="Freemind_Link_312753226" MODIFIED="1226225068760" TEXT="unique molecules but not necessarily biological molecules"/>
<node CREATED="1226225070857" ID="Freemind_Link_1322798370" MODIFIED="1226225072871" TEXT="maps to">
<node CREATED="1226225077010" ID="Freemind_Link_469758624" MODIFIED="1226225080209" TEXT="CAS"/>
<node CREATED="1226225080949" ID="Freemind_Link_1365093243" MODIFIED="1226225082869" TEXT="Smiles"/>
<node CREATED="1226225083408" ID="Freemind_Link_377956420" MODIFIED="1226225086562" TEXT="inChi"/>
<node CREATED="1226225087102" ID="Freemind_Link_1762459107" MODIFIED="1226225108632" TEXT="synonyms"/>
<node CREATED="1226225109171" ID="Freemind_Link_1774707598" MODIFIED="1226225111542" TEXT="platform"/>
<node CREATED="1226225120416" ID="Freemind_Link_1456086510" MODIFIED="1226225132263" TEXT="original name"/>
<node CREATED="1226225178187" ID="Freemind_Link_1179656018" MODIFIED="1226225179062" TEXT="RI"/>
</node>
</node>
<node CREATED="1226225008744" ID="Freemind_Link_340518378" MODIFIED="1226225012403" TEXT="BioMet">
<node CREATED="1226225237766" ID="Freemind_Link_1452454917" MODIFIED="1226225257261" TEXT="unique molecules that could be found in in vivo"/>
<node CREATED="1226225258815" ID="Freemind_Link_304065794" MODIFIED="1226225260561" TEXT="maps to">
<node CREATED="1226225262251" ID="Freemind_Link_722765624" MODIFIED="1226225264839" TEXT="KEGG"/>
<node CREATED="1226225265539" ID="Freemind_Link_337991948" MODIFIED="1226225272096" TEXT="*cyc"/>
<node CREATED="1226225276857" ID="Freemind_Link_1890422117" MODIFIED="1226225293748" TEXT="possible parent Analyte"/>
</node>
</node>
</node>
<node CREATED="1226228350058" ID="Freemind_Link_840668149" MODIFIED="1226228377154" POSITION="right" TEXT="implementation">
<node CREATED="1226228383529" ID="Freemind_Link_1310890551" MODIFIED="1226228389892" TEXT="SQLite">
<node CREATED="1226228392488" ID="Freemind_Link_1487938418" MODIFIED="1226228435134" TEXT="easy and portable, single file, could be easily be put somewhere on the net and then accessed from other machines"/>
</node>
<node CREATED="1226228356790" ID="Freemind_Link_31041804" MODIFIED="1226228459572" TEXT="create">
<node CREATED="1226228468649" ID="Freemind_Link_611350353" MODIFIED="1226228482586" TEXT="a perl scripts for every file type"/>
</node>
</node>
</node>
</map>
