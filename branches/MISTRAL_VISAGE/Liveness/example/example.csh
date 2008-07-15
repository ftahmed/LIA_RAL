echo "============================="
echo "  Live recording experience  "
echo "============================="

../bin/liveness --audio-feat ./data/DVD1_param_son_001.txt --video-feat ./data/DVD1_param_video_001.txt --energy ./data/DVD1_param_energie_trame_001.txt 

echo
echo "============================="
echo " Replay attack experience 1  "
echo "============================="
../bin/liveness --audio-feat ./data/DVD1_param_son_001.txt --video-feat ./data/DVD1_param_video_000.txt --energy ./data/DVD1_param_energie_trame_001.txt

echo
echo "============================="
echo " Replay attack experience 2  "
echo "============================="
../bin/liveness --audio-feat ./data/DVD1_param_son_001.txt --video-feat ./data/DVD1_param_video_002.txt --energy ./data/DVD1_param_energie_trame_001.txt
