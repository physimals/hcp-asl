import sys

from scripts import run_pipeline

if __name__ == "__main__":
    subid = "HCD0378150_V1_MR"

    cmd = f"""     
        --subid {subid} 
        --subdir /Users/thomaskirk/Data/{subid}
        --mbpcasl /Users/thomaskirk/Data/{subid}/unprocessed/mbPCASLhr/{subid}_mbPCASLhr_PA.nii.gz  
        --fmap_ap /Users/thomaskirk/Data/{subid}/unprocessed/mbPCASLhr/{subid}_PCASLhr_SpinEchoFieldMap_AP.nii.gz  
        --fmap_pa /Users/thomaskirk/Data/{subid}/unprocessed/mbPCASLhr/{subid}_PCASLhr_SpinEchoFieldMap_PA.nii.gz 
        --grads coeff_AS82_Prisma.grad
        --cores 2
        --regname MSMSulc 
    """

    cmd = " ".join(cmd.split())
    print(cmd)
    sys.argv[1:] = cmd.split()
    run_pipeline.main()
