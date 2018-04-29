RNASeq_Basics


-Into the working directory : Place Rdas.
- USer defined:
############################################################################################
######### USER DEFINED #####################################################################
############################################################################################
#Your files need to have a name that contains your comparisons "Olineo" "OlineoFTO" etc

        # Directory where the FPKMs ARE
        MyFPKMDirectory<-"./1_FPKMS/"

        #The Comparisons you want to make. The "name" should be contained AS IS in ALL the files
        #you want to compare
        Comparisons<-c("MEFs","Olineo","CBFA")

        #The different conditions. Can contain more than comparison. Will only look at the conditions
        #that contains the comparison within, without the Rep
        #Example: If you want to compare "Cancercell_condition1
        #(Cancercell_condition1_rep1,Cancercell_condition1_rep2, Cancercell_condition1_rep3)
        #and Cancercell_condition2 and Cancercell_condition3"etc etc
        #Your comparison should be "Cancercell" and in condition you should have
        #Cancercell_condition1, Cancercell_condition2, Cancercell_conditon3

        Conditions <- c("MEFs_CBFA_kno", "MEFs_CBFA_het",
                    "Olineo_CBFA_ctrol","Olineo_CBFA_siRNA")

                    
