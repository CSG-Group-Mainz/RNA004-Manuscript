### BASECALL TIME 
ml dorado
ml samtools

pod5_IVT="RNA004_peripheral_blood/2024_02_28_401_24_IVT"
pod5_Direct="RNA004_peripheral_blood/2024_02_28_401_24"

echo "Basecall 1: Blood IVT"
time dorado basecaller \
    /raid/chhewel_analysis/rna004_130bps_sup@v3.0.1 \
    $pod5_IVT \
    -r \
    --device "cuda:0,cuda:1,cuda:2,cuda:3" \
    > RNA004_blood_IVT_basecall.ubam
echo "FINISHED: Basecall 1: Blood IVT"
echo "########################"
echo "Basecall 2: Blood IVT"
time dorado basecaller \
    /raid/chhewel_analysis/rna004_130bps_sup@v3.0.1 \
    $pod5_IVT \
    --modified-bases m6A_DRACH \
    -r \
    --device "cuda:0,cuda:1,cuda:2,cuda:3" \
    > RNA004_blood_IVT_basecall.ubam
echo "FINISHED: Basecall 2: Blood IVT"
echo "########################"
echo "Basecall 3: Blood IVT"
time dorado basecaller \
    /raid/chhewel_analysis/rna004_130bps_sup@v3.0.1 \
    $pod5_IVT \
    --modified-bases m6A_DRACH \
    --estimate-poly-a \
    -r \
    --device "cuda:0,cuda:1,cuda:2,cuda:3" \
    > RNA004_blood_IVT_basecall.ubam
echo "FINISHED: Basecall 3: Blood IVT"
echo "########################"
echo "Basecall 1: Blood Normal"
time dorado basecaller \
    /raid/chhewel_analysis/rna004_130bps_sup@v3.0.1 \
    $pod5_Direct \
    -r \
    --device "cuda:0,cuda:1,cuda:2,cuda:3" \
    > RNA004_blood_Direct_basecall.ubam
echo "FINISHED: Basecall 1: Blood Normal"
echo "########################"
echo "Basecall 2: Blood Normal"
time dorado basecaller \
    /raid/chhewel_analysis/rna004_130bps_sup@v3.0.1 \
    $pod5_Direct \
    --modified-bases m6A_DRACH \
    -r \
    --device "cuda:0,cuda:1,cuda:2,cuda:3" \
    > RNA004_blood_Direct_basecall.ubam
echo "FINISHED: Basecall 2: Blood Normal"
echo "########################"
echo "Basecall 3: Blood Normal"
time dorado basecaller \
    /raid/chhewel_analysis/rna004_130bps_sup@v3.0.1 \
    $pod5_Direct \
    --modified-bases m6A_DRACH \
    --estimate-poly-a \
    -r \
    --device "cuda:0,cuda:1,cuda:2,cuda:3" \
    > RNA004_blood_Direct_basecall.ubam
echo "FINISHED: Basecall 3: Blood Normal"

