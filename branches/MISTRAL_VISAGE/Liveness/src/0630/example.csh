echo "============================="
echo "  Live recording experience  "
echo "============================="

./liveness --audio-feat ../Extract_Param_Bac/DVD1_param_son_001.txt --video-feat ../Extract_Param_Bac/DVD1_param_video_001.txt --energy ../Extract_Param_Bac/DVD1_param_energie_trame_001.txt 

echo
echo "============================="
echo " Replay attack experience 1  "
echo "============================="
./liveness --audio-feat ../Extract_Param_Bac/DVD1_param_son_001.txt --video-feat ../Extract_Param_Bac/DVD1_param_video_000.txt --energy ../Extract_Param_Bac/DVD1_param_energie_trame_001.txt

echo
echo "============================="
echo " Replay attack experience 2  "
echo "============================="
./liveness --audio-feat ../Extract_Param_Bac/DVD1_param_son_001.txt --video-feat ../Extract_Param_Bac/DVD1_param_video_002.txt --energy ../Extract_Param_Bac/DVD1_param_energie_trame_001.txt
