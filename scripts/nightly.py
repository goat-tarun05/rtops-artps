import sys, csv, os
sys.path.insert(0, ".")
from rtops.experiment import run_variant

def main(date):
    seed = int(date.replace("-", ""))
    r = run_variant({"session.date": date, "session.seed": seed}, out=f"runs/{date}")
    res, R, tcs = r["results"], r["reduced"], r["tcs"]
    
    pr1_after = res["pointing_runs"][0]["map_rms_true_after_deg"]
    pr2_resid = res["pointing_runs"][1]["rms_resid_before_deg"]
    sat_eff = R.get("sat", {}).get("mean_eff", 1.0)
    
    # Health Checks & Alarms
    alarms = []
    if pr2_resid > 0.05:
        alarms.append("PR2_High")
    if sat_eff < 0.7:
        alarms.append("SatEff_Low")
        
    for ob in tcs.obs:
        if ob.status == "SKIPPED":
            alarms.append(f"SKIPPED_{ob.id}")
            
    # Check last log status for each block to catch un-resumed wind aborts
    last_status = {seg["id"]: seg["status"] for seg in res["ob_log"]}
    for ob_id, status in last_status.items():
        if status == "ABORTED_WIND":
            alarms.append(f"WIND_{ob_id}")
            
    for t in res["tracking"]:
        if t.get("rms_true_deg", 0) > 0.3:
            alarms.append(f"TRK_{t.get('ob', 'Unknown')}")

    alarm_str = "|".join(alarms) if alarms else "OK"
    
    row = dict(date=date, seed=seed, 
               pr1_after=round(pr1_after, 4), 
               pr2_resid=round(pr2_resid, 4),
               sat_eff=round(sat_eff, 3) if "sat" in R else "",
               stows=len(res["stow_events"]), 
               skipped=sum(1 for o in tcs.obs if o.status == "SKIPPED"),
               alarms=alarm_str)
               
    os.makedirs("runs", exist_ok=True)
    new = not os.path.exists("runs/log.csv")
    with open("runs/log.csv", "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        if new: w.writeheader()
        w.writerow(row)
        
    print(f"[{date}] Completed. Alarms: {alarm_str}")
    if alarms:
        sys.exit(1)

if __name__ == "__main__":
    main(sys.argv[1])