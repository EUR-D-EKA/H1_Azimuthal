

int MergeTree(){
    TFile * input1 = TFile::Open("test_0.root","read");
    TFile * input2 = TFile::Open("test_1.root","read");
    TFile * input3 = TFile::Open("test_2.root","read");

    TTree *tree1 = (TTree*)input1->Get("properties");
    TTree *tree2 = (TTree*)input2->Get("properties");
    TTree *tree3 = (TTree*)input3->Get("properties"); //pointers to your 3 Trees


    TList *list = new TList;
    list->Add(tree1);
    list->Add(tree2);
    list->Add(tree3);
    TTree *newtree = TTree::MergeTrees(list);

    TFile *out = TFile::Open("out.root", "RECREATE");
    newtree->Write();
}