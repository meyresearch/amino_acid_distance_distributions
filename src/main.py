'''
Main 

'''
import functions
import PCN


if __name__ == '__main__':
    pcn = PCN.PCN('/Users/jasminguven/Downloads/2GIV.pdb')
    C_alphas = pcn.get_C_alphas()
   # pcn.link_length_threshold = 0
    ll = pcn.get_link_lengths(C_alphas)
    print(ll)
    #functions.main()