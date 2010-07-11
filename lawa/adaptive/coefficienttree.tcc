/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

namespace lawa {

// Is mu a descendant of lambda?
template <typename T, Construction Cons >
bool
isDescendant(const WaveletIndex<T,Cons> &lambda,
             const WaveletIndex<T,Cons> &mu)
{
    if ((lambda.j>mu.j) || 
        (lambda.j==mu.j && (lambda.xtype==mu.xtype || lambda.xtype==XWavelet))) {
        return false;
    }
	Support<T> BoxLambda=lambda.descendentCube(),
               BoxMu=mu.descendentCube();
    if (BoxLambda.l1<=BoxMu.l1 && BoxMu.l2<=BoxLambda.l2) {
        return true;
    } else {
        return false;
    }
}

template <typename T,Construction Cons >
Coefficient<Lexicographical,T,Cons>
completeTree(const Coefficient<Lexicographical,T,Cons> &coeff)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    const Basis<T,Primal,Interval,Cons> &basis=coeff.basis;
    Coefficient<Lexicographical,T,Cons> ret=coeff;

    for (it lambda=coeff.begin(); lambda!=coeff.end(); ++lambda) {
        if ((*lambda).first.xtype==XWavelet) {
            const WaveletIndex<T,Cons> &l=(*lambda).first;
            for (int j=l.j-1; j>=basis.j0-1; --j) {
                if (j==basis.j0-1) {
                    Support<T> BoxLambda=l.descendentCube();
                    for (WaveletIndex<T,Cons> mu(basis); mu.xtype==XBSpline; ++mu) {
                        Support<T> BoxMu=mu.descendentCube();
                        if (BoxMu.l1<=BoxLambda.l1 && BoxLambda.l2<=BoxMu.l2) {
                            if (coeff.count(mu)==0) {
                                ret.insert(val_type(mu,0.0));
                            }
                        }
                    }
                } else {
                    WaveletIndex<T,Cons> mu(basis,j,iceil(l.k*pow2i(-(l.j-j))),XWavelet);
                    if (coeff.count(mu)==0) {
                        ret.insert(val_type(mu,0.0));
                    }
                }
            }
        }
    }
    return ret;
}

template <typename T,Construction Cons >
tree<WaveletIndex<T,Cons> >
convert(const Coefficient<Lexicographical,T,Cons> &coeff)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator coeff_it;
    typedef typename tree<WaveletIndex<T,Cons> >::iterator tree_it;
    tree<WaveletIndex<T,Cons> > ret;
    for(coeff_it lambda=coeff.begin(); lambda!=coeff.end(); ++lambda){
        bool isSet=false;
        tree_it insert_node=ret.begin();
        for (tree_it leaf=ret.begin(); leaf!=ret.end(); ++leaf) {
            if (isDescendant(*leaf,(*lambda).first)) {
                if (ret.depth(insert_node)<=ret.depth(leaf)) {
                    insert_node=leaf;
                    isSet=true;
                }
            }
        }
        if (isSet==false) {
            ret.insert(ret.begin(), (*lambda).first);
        } else {
            ret.append_child(insert_node, (*lambda).first);
        }
    }
    return ret;
}

template <typename T,Construction Cons >
void
addToTree(tree<WaveletIndex<T,Cons> > &tr, const WaveletIndex<T,Cons> &index)
{
    typedef typename tree<WaveletIndex<T,Cons> >::iterator tree_it;
    bool isSet=false;
    tree_it insert_node=tr.begin();
    for (tree_it leaf=tr.begin(); leaf!=tr.end(); ++leaf) {
        if (isDescendant(*leaf,index)) {
            if (tr.depth(insert_node)<=tr.depth(leaf)) {
                insert_node=leaf;
                isSet=true;
            }
        }
    }
    if (isSet==false) {
        tr.insert(tr.begin(), index);
    } else {
        tr.append_child(insert_node, index);
    }
}

template <typename T,Construction Cons >
typename tree<WaveletIndex<T,Cons> >::iterator
findInTree(const tree<WaveletIndex<T,Cons> > &tr, const WaveletIndex<T,Cons> &index)
{
    typedef typename tree<WaveletIndex<T,Cons> >::fixed_depth_iterator fd_it;
    typedef typename tree<WaveletIndex<T,Cons> >::sibling_iterator sibl_it;
    typedef typename tree<WaveletIndex<T,Cons> >::iterator tree_it;
    // first, check all root nodes...
    tree_it it;
    for (fd_it root=tr.begin_fixed(tr.begin(),0); root!=tr.end_fixed(tr.begin(),0); ++root) {
        if (isDescendant(*root,index)) {
            it=root;
            break;
        } else if (*root==index) {
            return root;
        }
    }
    // ...then traverse through tree and search index
    do {
        for (sibl_it child=tr.begin(it); child!=tr.end(it); ++child) {
            if (isDescendant(*child,index)) {
                it=child;
                break;
            } else if (*child==index) {
                return child;
            }
        }
    } while (*it!=index && it!=tr.end());
    assert(*it==index);
    return it;
}

template <typename T,Construction Cons>
std::ostream&
operator<<(std::ostream &s, const tree<WaveletIndex<T,Cons> > &_tr)
{
    typedef typename tree<WaveletIndex<T,Cons> >::iterator tree_it;
    for (tree_it it=_tr.begin(); it!=_tr.end(); ++it) {
        for (int i=0; i<_tr.depth(it); ++i) {
		    s << "  ";
        }
        s << (*it) << endl;
    }
    return s;
}

} // namespace lawa
