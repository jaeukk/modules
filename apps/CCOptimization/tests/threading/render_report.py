"""Render the final Markdown report and interactive measurement chart as standalone HTML.

Requires pandoc. Run from any directory; inputs are relative to this script.
"""
import html
import json
from pathlib import Path
import subprocess


def main():
    app = Path(__file__).resolve().parents[2]
    records = {}
    for folder in sorted((Path(__file__).parent / "results-20261019").iterdir()):
        if folder.is_dir():
            records[folder.name] = {
                path.stem: json.loads(path.read_text()) for path in sorted(folder.glob("*.json"))
            }
    body = subprocess.check_output(
        ["pandoc", "-f", "gfm", "-t", "html", str(app / "THREADING_RESULTS.md")], text=True
    )
    title = "CCOptimization: scheduling and L-BFGS history"
    document = """<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>TITLE</title><style>
*{box-sizing:border-box}body{margin:0;background:#f5f7fb;color:#18283d;font:17px/1.65 system-ui,sans-serif}
main{max-width:1040px;margin:auto;padding:28px 24px}h1,h2{line-height:1.25}h2{margin-top:2.2em}
a{color:#175aa2}table{width:100%;border-collapse:collapse;font-size:.9em;display:block;overflow:auto}
td,th{padding:9px;border-bottom:1px solid #ced7e3;text-align:left;white-space:nowrap}
pre{overflow:auto;padding:16px;background:#e6ecf5}code{font-size:.87em}section.chart{padding:22px;background:white;border:1px solid #ccd7e8;border-radius:12px}
select,button{font:inherit;max-width:100%;padding:7px;margin:4px 8px 10px 0}
.pair{margin:15px 0}.track{display:flex;align-items:center;gap:8px;margin:4px 0}.bar{min-width:2px;height:18px;border-radius:3px;background:#2869b2}.b{background:#c46b12}
.value{white-space:nowrap;font-size:.85em}.key{font-weight:700;color:#2869b2}.key.b{background:none;color:#a45706}
@media(max-width:600px){main{padding:18px 14px}body{font-size:15px}section.chart{padding:12px}.value{font-size:.8em}}
</style><main><section class="chart" aria-labelledby="chart-title">
<h1 id="chart-title">Measured scheduling performance</h1>
<p>Three species · four bound CPUs · contact strength <strong>100</strong> · history 100.<br>
<span class="key">A: concurrent species</span> · <span class="key b">B: sequential species, internally threaded FINUFFT</span></p>
<label for="metric">Measurement</label><br><select id="metric">
<option value="eval">Warm evaluation wall time (ms/call)</option>
<option value="opt">Bounded NLopt wall time (seconds)</option>
<option value="cpu">Warm evaluation CPU utilization (%)</option></select>
<button id="download">Download measurement JSON</button><div id="bars" aria-live="polite"></div>
<p><strong>Zero accepted ground states in the A/B timing runs.</strong> These timings compare bounded workloads.
The accepted-ground-state history test below used contact strength 1; production stays at 100.</p></section>
BODY</main><script id="measurements" type="application/json">DATA</script><script>
'use strict';
const data=JSON.parse(document.getElementById('measurements').textContent);
const measurements=[...data['ab-small'].summary,...data['ab-4000'].summary];
function draw(){
 const metric=document.getElementById('metric').value;
 const rows=measurements.filter(x=>x.mode===(metric==='opt'?'opt':'eval'));
 const value=x=>metric==='cpu'?x.utilization_percent:metric==='eval'?x.wall_seconds*1000/400:x.wall_seconds;
 const maximum=metric==='cpu'?100:Math.max(...rows.map(value));
 const container=document.getElementById('bars');container.replaceChildren();
 for(const size of [300,600,1000,4000]){
  const pair=document.createElement('div');pair.className='pair';
  const label=document.createElement('strong');label.textContent=size.toLocaleString()+' particles/species';pair.append(label);
  for(const option of ['A','B']){
   const row=rows.find(x=>x.particles_per_species===size&&x.option===option);
   const track=document.createElement('div');track.className='track';
   const bar=document.createElement('div');bar.className='bar'+(option==='B'?' b':'');bar.style.width=(65*value(row)/maximum)+'%';
   const text=document.createElement('span');text.className='value';text.textContent=option+': '+value(row).toFixed(metric==='cpu'?1:3)+(metric==='cpu'?'%':metric==='eval'?' ms/call':' s');
   track.append(bar,text);pair.append(track);
  }container.append(pair);
 }
}
document.getElementById('metric').addEventListener('change',draw);
document.getElementById('download').addEventListener('click',()=>{
 const url=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)],{type:'application/json'}));
 const link=document.createElement('a');link.href=url;link.download='cco-threading-history-20261019.json';link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
});draw();
</script></html>"""
    document = document.replace("TITLE", html.escape(title)).replace("BODY", body)
    document = document.replace("DATA", json.dumps(records).replace("<", "\\u003c"))
    (app / "THREADING_RESULTS.html").write_text(document)


if __name__ == "__main__":
    main()
