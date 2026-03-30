#!/usr/bin/expect -f

set timeout -1

# Launch HyPhy interactively
spawn hyphy

# Send input matching HyPhy's prompts
expect "Please select type of analyses*" {
    send "13\r"
}

expect "Please select the analysis*" {
    send "3\r"
}

expect "Please choose an option*" {
    send "1\r"
}

expect "Please specify a codon data file*" {
    send "/scratch/sguirale/H5N1_Antivirals/PB1/raxml2/PB1_aln_pruned_mesquite.phy\r"
}

expect "Please select a tree file*" {
    send "/scratch/sguirale/H5N1_Antivirals/PB1/raxml2/RAxML_bestTree.PB1.tree\r"
}

expect "Please choose an option*" {
    send "2\r"
}

expect "Please choose option*" {
    send "9\r"
}

expect "Please choose option*" {
    send "d\r"
}

expect "Please choose an option*" {
    send "1\r"
}

expect "Please choose an option*" {
    send "1\r"
}

expect "Choose the cutoff*" {
    send "0.95\r"
}

expect "Choose the number of categories*" {
    send "2\r"
}

expect "Write detailed results*" {
    send "/scratch/sguirale/H5N1_Antivirals/PB1/hyphy/PB1.resultsM8.json\r"
}

# Wait for completion
expect eof
