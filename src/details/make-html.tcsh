#!/bin/tcsh -ef
# author : sjn
# date : Fri Sep 23 12:29:05 PDT 2022

# This script assumes the qc graphics are in the same dir as the <html-out>

if ( $#argv < 4 ) then
  printf "Expect $0 <sample-name> <html-overview> <html-meat> <pdf|txt>*\n"
  exit -1
endif

set echo

set sample=$1
set overview=$2 # also an html output
set html=$3
set pdfs = "$argv[4-]"
set txts = `echo "$pdfs" | tr ' ' '\n' | grep -e ".txt"`
set pdfs = `echo "$pdfs" | tr ' ' '\n' | grep -e ".pdf"`

cat <<__HTML__ >! $html
<html>
  <head>
    <style>
      table, th, td {
        border:1px solid black;
      }
      img {
        border: 1px solid #ddd;
        border-radius: 4px;
        padding: 5px;
        width: 150px;
      }
      img:hover {
        box-shadow: 0 0 2px 1px rgba(0, 140, 186, 0.5);
      }
    </style>
  </head>
  <body>
    <table>
    <caption><h1>QC Report: $sample</h1></caption>
__HTML__

@ cntr = 1
foreach pdf ($pdfs)
  convert $pdf $pdf:r.png
  set nm = `echo $pdf:t:r | sed 's;qc_;;g'`
  printf '      <tr>' >>! $html
  printf '        <th>'$nm'</th>' >>! $html
  printf '        <td style="text-align: center;"><img src="'$pdf:t:r.png'" alt="'$pdf:t'" onmouseenter="hover(this)"></img></td>' >>! $html
  printf '        <td style="text-align: center;"><pre>' >>! $html
  cat $txts[$cntr] | grep -vi "# Note: \*" >>! $html
  printf '        </pre></td>\n' >>! $html
  printf '      </tr>\n' >>! $html
  @ cntr++
end

cat << __HTML__ >>! $html
    </table>
    <script>
      function hover(element) {
        window.parent.postMessage(element.alt, '*')
      }
    </script>
  </body>
</html>
__HTML__


set pdf = `echo $pdfs | tr ' ' '\n' | grep -e ".pdf" | head -1`
cat <<__HTML__ >! $overview
<html>
  <head>
    <style>
      html,body {
        width: 100%;
        height: 100%;
      }
    </style>
  </head>
  <body>
    <div class="container">
      <iframe id="TNAILS" src="$html:t" frameborder="0" scrolling="yes" style="height: 100%; width: 45%; float: left;" align="left"></iframe>
      <iframe id="FULLIMG" src="$pdf:t" frameborder="0" scrolling="yes" style="overflow: hidden; height: 100%; width:55%;" align="right"></iframe>
      <script>
        window.onmessage = function(msg) {
          const img = document.getElementById("FULLIMG")
          img.setAttribute('src', msg.data)
        };
      </script>
    </div>
  </body>
</html>
__HTML__

exit 0
