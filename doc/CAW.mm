<map version="0.9.0">
<!-- To view this file, download free mind mapping software FreeMind from http://freemind.sourceforge.net -->
<node CREATED="1457342480282" ID="ID_1224971888" MODIFIED="1457431980543" TEXT="CAW">
<node CREATED="1457346470715" FOLDED="true" ID="ID_1659796874" MODIFIED="1458031779282" POSITION="right" TEXT="Prototyping">
<node CREATED="1457346528430" FOLDED="true" ID="ID_1279764966" MODIFIED="1457436193849" TEXT="Pall&apos;s make">
<node CREATED="1457347260782" ID="ID_1444683302" MODIFIED="1457347634840">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      In this prototype we already have
    </p>
    <ul>
      <li>
        <b>gmake</b>&#160;as a general language used by all UNIX
      </li>
      <li>
        GATK best practices implemented that can run on a single huge machine
      </li>
      <li>
        using a restricted set of <b>make</b>&#160;to prevent code from clutter
      </li>
      <li>
        easy to use for a single person
      </li>
    </ul>
    <p>
      Missing parts
    </p>
    <ul>
      <li>
        not optimized for cluster/LIMS use
      </li>
      <li>
        can be tricky to keep is running for hundreds of samples&#160;
      </li>
    </ul>
  </body>
</html></richcontent>
</node>
</node>
<node CREATED="1457346545122" ID="ID_315494008" MODIFIED="1457527404988" TEXT="Malin&apos;s shell">
<node CREATED="1457347650139" ID="ID_1348560814" MODIFIED="1457427579784">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      In this prototype we already have
    </p>
    <ul>
      <li>
        T/N workflow implemented in bash (also general UNIX language)
      </li>
      <li>
        additional functionalities added that are not included in GATK BP
      </li>
      <li>
        Version control set up at Bitbucket
      </li>
      <li>
        slurm included
      </li>
      <li>
        a comprehensive prototype showing many needs by practicing researchers
      </li>
    </ul>
    <p>
      Missing parts
    </p>
    <ul>
      <li>
        do not know frankly, waiting for input from Malin
      </li>
      <li>
        thorough testing and db connection maybe
      </li>
    </ul>
  </body>
</html></richcontent>
</node>
</node>
<node CREATED="1457346555370" ID="ID_761801786" MODIFIED="1457434805298" TEXT="other canned things on the net">
<node CREATED="1457348602389" FOLDED="true" ID="ID_407753552" MODIFIED="1457436197497" TEXT="BCBIO">
<node CREATED="1457348646219" ID="ID_202427485" MODIFIED="1457348654855" TEXT="http://bcbio-nextgen.readthedocs.org/en/latest/contents/testing.html"/>
<node CREATED="1457348657358" ID="ID_193810972" MODIFIED="1457426013986">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      They have a pretty comprehensive working platform
    </p>
    <ul>
      <li>
        multiple VCs
      </li>
      <li>
        T/N comparison
      </li>
      <li>
        parallelisation IPython parallel
      </li>
      <li>
        ...
      </li>
    </ul>
    <p>
      Seems they went through the pains at least some level and they are usually pretty helpful
    </p>
  </body>
</html></richcontent>
</node>
</node>
<node CREATED="1457426525556" ID="ID_189905720" MODIFIED="1457426530910" TEXT="GATK CGA"/>
<node CREATED="1457426650352" ID="ID_733580872" MODIFIED="1457427783714" TEXT="and many others to be reviewed ..."/>
</node>
</node>
<node CREATED="1457346489629" FOLDED="true" ID="ID_1905157975" MODIFIED="1458031776927" POSITION="left" TEXT="Final Workflow">
<node CREATED="1457429715825" ID="ID_1887772963" MODIFIED="1457524852116" TEXT="Agreed to be written in nextflow (Szilva)">
<icon BUILTIN="full-1"/>
</node>
<node CREATED="1457429742188" ID="ID_1621399207" MODIFIED="1457524699195" TEXT="GenStats integration">
<icon BUILTIN="full-2"/>
</node>
<node CREATED="1457429729459" ID="ID_397088719" MODIFIED="1457524838653" TEXT="include into NGI-PL (Szilva)">
<icon BUILTIN="full-2"/>
</node>
<node CREATED="1457430981253" ID="ID_1107801858" MODIFIED="1457430986643" TEXT="Integrate MQC"/>
<node CREATED="1457431100904" ID="ID_1056710847" MODIFIED="1457431910057" TEXT="Things to deliver">
<node CREATED="1457431151357" ID="ID_766837473" MODIFIED="1457436303233" TEXT="VCF">
<icon BUILTIN="full-1"/>
</node>
<node CREATED="1457431159942" ID="ID_1234908553" MODIFIED="1457436307538" TEXT="alignments">
<icon BUILTIN="full-1"/>
<node CREATED="1457437916222" ID="ID_178827813" MODIFIED="1457437932620" TEXT="bwa mem parameters (Pall &amp; Malin)"/>
<node CREATED="1457437937754" ID="ID_51402994" MODIFIED="1457437976789" TEXT="base realignment &amp; recal &amp; MD (Pall &amp; Malin)"/>
</node>
<node CREATED="1457431165638" ID="ID_428674581" MODIFIED="1457436311850" TEXT="QC">
<icon BUILTIN="full-1"/>
</node>
<node CREATED="1457431911187" ID="ID_391994940" MODIFIED="1457527220079" TEXT="Annotation to be postponed for now - use Gemini">
<icon BUILTIN="full-2"/>
</node>
<node CREATED="1457437373096" ID="ID_1202295689" MODIFIED="1457524770285" TEXT="Structural Variants (Pelin)">
<icon BUILTIN="full-2"/>
</node>
</node>
<node CREATED="1457431109936" FOLDED="true" ID="ID_284429603" MODIFIED="1457524690212" TEXT="Things NOT to deliver">
<node CREATED="1457431191215" ID="ID_185666443" MODIFIED="1457431272134" TEXT="Clinically actionable stuff"/>
</node>
<node CREATED="1457436798044" ID="ID_21929598" MODIFIED="1457525687648" TEXT="Tasks to do">
<icon BUILTIN="xmag"/>
<node CREATED="1457436821600" ID="ID_668437699" MODIFIED="1457524804814" TEXT="Install VarDict on UPPMAX">
<icon BUILTIN="full-2"/>
</node>
<node CREATED="1457436832770" ID="ID_1856808140" MODIFIED="1457524687106" TEXT="Look at MuTect2 speed issues (Szilva)">
<icon BUILTIN="full-2"/>
</node>
<node CREATED="1457436847286" ID="ID_409574782" MODIFIED="1457524681683" TEXT="benchmarking for disk usage">
<icon BUILTIN="full-3"/>
</node>
<node CREATED="1457437115254" ID="ID_108247371" MODIFIED="1457527295415" TEXT="distributing parallel VC processes sctratch per node (run in one node or how to distribute)">
<icon BUILTIN="full-2"/>
</node>
<node CREATED="1457437215832" FOLDED="true" ID="ID_543463510" MODIFIED="1457527292055" TEXT="SV &amp; CNV call (Pelin)">
<icon BUILTIN="full-1"/>
<node CREATED="1457437240431" ID="ID_1760728208" MODIFIED="1457437257486" TEXT="meeting to Francesco and Jesper"/>
<node CREATED="1457437260036" ID="ID_1592833371" MODIFIED="1457437280940" TEXT="presenting their stuff"/>
</node>
<node CREATED="1457438358320" ID="ID_780924236" MODIFIED="1457527280423" TEXT="Merging varinats (SNP &amp; indels only)">
<icon BUILTIN="full-1"/>
</node>
</node>
</node>
<node CREATED="1457346508046" FOLDED="true" ID="ID_568515110" MODIFIED="1458031780612" POSITION="right" TEXT="Version control">
<icon BUILTIN="full-1"/>
<node CREATED="1457426676542" ID="ID_1303068604" MODIFIED="1457426687774" TEXT="Seems we settled with Github"/>
<node CREATED="1457426690481" ID="ID_472426782" MODIFIED="1457435269899">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      Have to consider that things we are treating as private now witll be eventually public
    </p>
    <ul>
      <li>
        private repos are available for either SciLifeLab or NGI for free
      </li>
      <li>
        deployment should happen from Git to allow continuous integration&#160;
      </li>
      <li>
        ask Malin how she feels about being the responsible person for Github
      </li>
    </ul>
  </body>
</html></richcontent>
</node>
<node CREATED="1457430450034" ID="ID_960477707" MODIFIED="1457526098181" TEXT="Issue tracking: use Trello"/>
</node>
<node CREATED="1457346603206" FOLDED="true" ID="ID_106204864" MODIFIED="1457528238652" POSITION="left" TEXT="Documentation">
<icon BUILTIN="button_ok"/>
<node CREATED="1457430380443" ID="ID_287688980" MODIFIED="1457430395910" TEXT="Will need a comprehensive user manual"/>
<node CREATED="1457430405721" ID="ID_1444277685" MODIFIED="1457430417303" TEXT="Articles related to methods used">
<node CREATED="1457432019871" ID="ID_1153709687" MODIFIED="1457435145138" TEXT="http://www.nature.com/ncomms/2015/151209/ncomms10001/abs/ncomms10001.html"/>
<node CREATED="1457432116935" ID="ID_1215594293" MODIFIED="1457432154934" TEXT="http://www.ncbi.nlm.nih.gov/pubmed/25984700"/>
</node>
<node CREATED="1457430421914" ID="ID_316876813" MODIFIED="1457430432676" TEXT="Everything should be stored on Github"/>
<node CREATED="1457431337292" ID="ID_1840623386" MODIFIED="1457431418586" TEXT="https://wabi-wiki.scilifelab.se/display/SHGATG/Somatic+variant+calling+in+cancer"/>
</node>
<node CREATED="1457346617090" ID="ID_1466099917" MODIFIED="1458031745463" POSITION="right" TEXT="Sample data">
<node CREATED="1457427788368" ID="ID_370049090" MODIFIED="1457427796216" TEXT="DREAM challenge"/>
<node CREATED="1457427805718" ID="ID_568969957" MODIFIED="1457427814118" TEXT="for workflow testing"/>
<node CREATED="1457427816067" ID="ID_296232360" MODIFIED="1457524554235" TEXT="Benchmark article before Xmas with public data? (somebody was mentioning, but which article is it?)"/>
<node CREATED="1457429648581" ID="ID_842175135" MODIFIED="1457429658021" TEXT="Do we have real validated data?"/>
</node>
<node CREATED="1457346622414" FOLDED="true" ID="ID_1693603794" MODIFIED="1457524588671" POSITION="left" TEXT="Real Data">
<icon BUILTIN="button_ok"/>
<node CREATED="1457429666520" ID="ID_1951908802" MODIFIED="1457437776669" TEXT="Minor issues">
<node CREATED="1457430244990" ID="ID_1641102756" MODIFIED="1457430247495" TEXT="XML2yaml in NGI-PL"/>
<node CREATED="1457430249688" ID="ID_1611273069" MODIFIED="1457430275302" TEXT="Pulling IDs from Charon"/>
</node>
<node CREATED="1457430124931" ID="ID_1502450940" MODIFIED="1457430149959" TEXT="Some real data will be needed that we can use for development"/>
<node CREATED="1457430152385" ID="ID_236813007" MODIFIED="1457437786572" TEXT="What will we get? FASTQs only, or Illumina dir structure? (Szilva to ask Francesco)">
<icon BUILTIN="full-1"/>
</node>
<node CREATED="1457430288362" ID="ID_964932531" MODIFIED="1457434109831" TEXT="Data inflow">
<node CREATED="1457430309221" ID="ID_1375635014" MODIFIED="1457430315139" TEXT="90 T/N pairs first"/>
<node CREATED="1457430316842" ID="ID_722589" MODIFIED="1457430326855" TEXT="half a dozen already sequenced"/>
<node CREATED="1457430329269" ID="ID_1771384866" MODIFIED="1457434155513" TEXT="WES or WGS? It is expected to be WGS"/>
<node CREATED="1457430338844" ID="ID_1725851591" MODIFIED="1457430365443" TEXT="more data to come distributed between different projects in bursts"/>
</node>
<node CREATED="1457431020826" ID="ID_1961862583" MODIFIED="1457437795590" TEXT="QC">
<icon BUILTIN="button_ok"/>
<node CREATED="1457431026038" ID="ID_1819752136" MODIFIED="1457434361326" TEXT="MultiQC for bioinfo parts (Szilva &amp; Phil)">
<icon BUILTIN="full-3"/>
</node>
<node CREATED="1457431042487" ID="ID_1468332204" MODIFIED="1457434400750" TEXT="Sample heterogenity QC (Malin)">
<icon BUILTIN="full-3"/>
<node CREATED="1457431591080" ID="ID_1464729466" MODIFIED="1457431607232" TEXT="ASCAT in PNAS"/>
<node CREATED="1457431998154" ID="ID_1900280895" MODIFIED="1457432005245" TEXT="MATH"/>
</node>
<node CREATED="1457432199197" ID="ID_1756607997" MODIFIED="1457434763328" TEXT="Sequence specific sequencing errors (sample prep related)"/>
<node CREATED="1457434547693" ID="ID_296241543" MODIFIED="1457434573551" TEXT="QC3 (Szilva)">
<icon BUILTIN="full-3"/>
</node>
</node>
</node>
<node CREATED="1457347056851" ID="ID_883787027" MODIFIED="1458031743923" POSITION="right" TEXT="Validation">
<node CREATED="1457429766837" ID="ID_886309607" MODIFIED="1457429786012" TEXT="Mostly the VC that needs to be validated"/>
<node CREATED="1457429912608" ID="ID_1519835453" MODIFIED="1457429957618" TEXT="Look at the alternatives used by other groups"/>
<node CREATED="1457429842381" ID="ID_706834882" MODIFIED="1457429849043" TEXT="VCs to be included">
<node CREATED="1457429850540" ID="ID_49374550" MODIFIED="1457429857393" TEXT="MuTect2"/>
<node CREATED="1457429859605" ID="ID_1964895612" MODIFIED="1457429862833" TEXT="FreeBayes"/>
<node CREATED="1457429865029" ID="ID_534792796" MODIFIED="1457429940830" TEXT="SV Jesper ?"/>
</node>
<node CREATED="1457430066261" ID="ID_930004847" MODIFIED="1457430093479" TEXT="Do we have criteria for QC here?"/>
<node CREATED="1457430096638" ID="ID_1685082614" MODIFIED="1457430110502" TEXT="Do we have real data to validate to?"/>
</node>
<node CREATED="1457430718213" ID="ID_326519709" MODIFIED="1458034084893" POSITION="left" TEXT="People involved">
<icon BUILTIN="button_ok"/>
<node CREATED="1457430909753" ID="ID_1312483697" MODIFIED="1457430915044" TEXT="SciLifeLab">
<node CREATED="1457430920150" ID="ID_1509192982" MODIFIED="1457438084705" TEXT="Valtteri Wirta (Head of Clin. Seq)"/>
<node CREATED="1457430726521" ID="ID_686345913" MODIFIED="1457438148366" TEXT="Pelin ">
<icon BUILTIN="button_ok"/>
</node>
<node CREATED="1457430763142" ID="ID_1341216236" MODIFIED="1457436132409" TEXT="Max"/>
<node CREATED="1457430767794" ID="ID_875113290" MODIFIED="1457438151478" TEXT="Szilva">
<icon BUILTIN="button_ok"/>
</node>
<node CREATED="1457430772162" ID="ID_816704477" MODIFIED="1457438155750" TEXT="Pall">
<icon BUILTIN="button_ok"/>
</node>
<node CREATED="1457430777643" ID="ID_1055740395" MODIFIED="1457438159166" TEXT="Malin">
<icon BUILTIN="button_ok"/>
</node>
<node CREATED="1457430782434" ID="ID_904080224" MODIFIED="1457436119632" TEXT="Bj&#xf6;rn"/>
<node CREATED="1457430802877" ID="ID_1546594007" MODIFIED="1457436123312" TEXT="Phil (as a consultant)"/>
<node CREATED="1457433920538" ID="ID_957795875" MODIFIED="1458034099215" TEXT="P&#xe4;r Lundin (as a consultant)">
<icon BUILTIN="help"/>
</node>
<node CREATED="1457434956438" ID="ID_1040934830" MODIFIED="1457434972415" TEXT="Francesco (as a consultant)"/>
<node CREATED="1457437617848" ID="ID_371820312" MODIFIED="1457438164742" TEXT="Jesper">
<icon BUILTIN="button_ok"/>
</node>
</node>
<node CREATED="1457430872537" ID="ID_1316963676" MODIFIED="1457430875163" TEXT="KI">
<node CREATED="1457430876254" ID="ID_1041384326" MODIFIED="1457436145937" TEXT="Monica Nister"/>
<node CREATED="1457430880855" ID="ID_926709286" MODIFIED="1457436142329" TEXT="Teresita D&#xed;az De St&#xe5;hl"/>
</node>
</node>
<node CREATED="1457431534390" ID="ID_1265262821" MODIFIED="1458031742322" POSITION="right" TEXT="Reference">
<node CREATED="1457431543875" ID="ID_1723084016" MODIFIED="1457431547612" TEXT="HG19"/>
<node CREATED="1457431553203" ID="ID_5179918" MODIFIED="1457431569987" TEXT="GRCh38 aka HG20"/>
<node CREATED="1457431629851" ID="ID_1066979485" MODIFIED="1457431639290" TEXT="HG19+decoy"/>
</node>
<node CREATED="1457525760313" ID="ID_1787073667" MODIFIED="1458034220909" POSITION="left" TEXT="Spring 1">
<node CREATED="1457525779915" ID="ID_1640138259" MODIFIED="1457526135327" TEXT="week 1">
<node CREATED="1457525802881" ID="ID_1075735433" MODIFIED="1458034333803" TEXT="Set up management tools (Szilva)">
<node CREATED="1457525924169" ID="ID_1195720276" MODIFIED="1457525936283" TEXT="Github"/>
<node CREATED="1457525939813" ID="ID_387404673" MODIFIED="1457525942813" TEXT="Trello"/>
<node CREATED="1458034335242" ID="ID_1031306725" MODIFIED="1458034350615" TEXT="Include project plan in Google doc">
<node CREATED="1458034351997" ID="ID_596766653" MODIFIED="1458034351997" TEXT=""/>
</node>
</node>
<node CREATED="1457525818905" ID="ID_1812505116" MODIFIED="1457600266903" TEXT="Finalize individual steps - have complete graph view (Malin,Pelin, Pall)">
<node CREATED="1457527427077" ID="ID_361455930" MODIFIED="1457527444521" TEXT="Small set of VSs"/>
<node CREATED="1457527447086" ID="ID_299668514" MODIFIED="1457527459083" TEXT="Merging VCFs"/>
<node CREATED="1457527469876" ID="ID_1583992897" MODIFIED="1457527485882" TEXT="alignment and recal steps"/>
</node>
<node CREATED="1457525836626" ID="ID_429492813" MODIFIED="1457526223284" TEXT="Sample DB, sample file structure (Szilva)"/>
<node CREATED="1457526166824" ID="ID_1095703251" MODIFIED="1457527060986" TEXT="Propose framework for SV (Pelin)"/>
<node CREATED="1457527369329" ID="ID_1656558706" MODIFIED="1457527725690" TEXT="propose framework for sample heterogenity QC (Malin)"/>
<node CREATED="1457528173301" ID="ID_1489526038" MODIFIED="1457528202278" TEXT="Learning NexFlow (Malin, Pall, Szilva)"/>
</node>
<node CREATED="1457525791476" ID="ID_227619460" MODIFIED="1458031785713" TEXT="week 2">
<node CREATED="1457527121458" ID="ID_1884939564" MODIFIED="1457527133792" TEXT="Set up test data">
<node CREATED="1457527143848" ID="ID_575526" MODIFIED="1457527166139" TEXT="for quick runs (i.e. unit testing)"/>
<node CREATED="1457527168929" ID="ID_420957985" MODIFIED="1457527196333" TEXT="set up data for system testing (DREAM?)"/>
</node>
<node CREATED="1457527318111" ID="ID_1793324672" MODIFIED="1457527557567" TEXT="propose initial QC"/>
<node CREATED="1457527611193" ID="ID_1566426139" MODIFIED="1457527628241" TEXT="test sample QC"/>
<node CREATED="1457527638215" ID="ID_554304240" MODIFIED="1457527641676" TEXT="test SV"/>
<node CREATED="1457527646168" ID="ID_1265704748" MODIFIED="1457527688873" TEXT="learn Github, Trello"/>
<node CREATED="1457527560421" ID="ID_762906328" MODIFIED="1458031933760" TEXT="implement basic steps in NF (Pall)">
<icon BUILTIN="button_ok"/>
<node CREATED="1458031956744" ID="ID_749031960" MODIFIED="1458031987307" TEXT="It is doing something, but not perfect yet"/>
</node>
</node>
<node CREATED="1457525909772" ID="ID_934800376" MODIFIED="1457525913888" TEXT="week 3"/>
</node>
</node>
</map>
