# 2sub.sh revision

ml python
python scripts/rep_frames.py -s nw -n new_clusters
python scripts/rep_frames.py -s nw -n new_clusters -t 0.1
python scripts/rep_frames.py -s nw -n new_clusters -t 0.2
python scripts/rep_frames.py -s nw -n new_clusters -t 0.3
python scripts/rep_frames.py -s nw -n new_clusters -i SM
python scripts/rep_frames.py -s nw -n new_clusters -t 0.1 -i SM
python scripts/rep_frames.py -s nw -n new_clusters -t 0.2 -i SM
python scripts/rep_frames.py -s nw -n new_clusters -t 0.3 -i SM

ml python
python scripts/rep_frames.py -s nw -n v3_norm
python scripts/rep_frames.py -s nw -n v3_norm -t 0.1
python scripts/rep_frames.py -s nw -n v3_norm -t 0.2
python scripts/rep_frames.py -s nw -n v3_norm -t 0.3
python scripts/rep_frames.py -s nw -n v3_norm -i SM
python scripts/rep_frames.py -s nw -n v3_norm -t 0.1 -i SM
python scripts/rep_frames.py -s nw -n v3_norm -t 0.2 -i SM
python scripts/rep_frames.py -s nw -n v3_norm -t 0.3 -i SM

#normalize sub.sh revision
python normalize.py -b 750 -n v3