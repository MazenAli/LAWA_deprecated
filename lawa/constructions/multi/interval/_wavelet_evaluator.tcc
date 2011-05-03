namespace lawa {

//--- linear evaluators --------------------------------------------------------

template <typename T>
T
_linear_wavelet_left_evaluator0(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (0.125 <= x && x < 0.25) {
            value = 8.42740770417694358377 - 40.7411946788691781932*x;
        } else if(0.5 <= x && x < 0.75){
            value = 1.174108337530388360919 - 1.643751672542543705286*x;
        } else if(0.75 <= x && x < 1.0){
            value = -0.2348216675060776721838 + 0.2348216675060776721838*x;
        } else if(0.0 <= x && x < 0.125){
            value = 26.67806695454637047700 * x;
        } else if(0.375 <= x && x <0.5){
            value = -2.92872776231550774859 + 6.56192052714924851372*x;
        } else if(0.25 <= x && x < 0.375){
            value = -4.33765776735197378169 + 10.31906720724649126866*x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0.125 <= x && x < 0.25) {
            value =  - 40.7411946788691781932;
        } else if(0.5 <= x && x < 0.75){
            value =  - 1.643751672542543705286;
        } else if(0.75 <= x && x < 1.0){
            value =  0.2348216675060776721838;
        } else if(0.0 <= x && x < 0.125){
            value = 26.67806695454637047700;
        } else if(0.375 <= x && x <0.5){
            value =  6.56192052714924851372;
        } else if(0.25 <= x && x < 0.375){
            value =  10.31906720724649126866;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_left_evaluator1(T x, unsigned short deriv)
{
    double value = 0.0;
    //x = power2(-1.0)*x;//anpasssung
    if (deriv == 0) {
        if (0.0 <= x && x < 0.125) {
            value = 5.2638582743553581630 - 53.356133909092740954*x;
        } else if (0.25 <= x && x < 0.375){
            value = 1.1741083375303883609 - 3.287503345085087411*x;
        } else if (0.375 <= x && x < 0.5){
            value = -0.2348216675060776722 + 0.469643335012155344*x;
        } else if (0.125 <= x && x < 0.25){
            value = -3.163549429821585421 + 14.063127724322807716*x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0.0 <= x && x < 0.125) {
            value = - 53.356133909092740954;
        } else if (0.25 <= x && x < 0.375){
            value =  - 3.287503345085087411;
        } else if (0.375 <= x && x < 0.5){
            value =  0.469643335012155344;
        } else if (0.125 <= x && x < 0.25){
            value =  14.063127724322807716;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_inner_evaluator0(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (0 <= x && x < 0.125) {
            value = 2.17124059336723766167 - 33.0125730435099442245 * x;
        } else if(-0.375 <= x && x < -0.25) {
            value = -7.14234036204314274860 - 19.69261569230427143433 * x;
        } else if(0.25 <= x && x < 0.375){
            value = 2.27349781934085223456 - 5.6124494201652150829 * x;
        } else if(-0.5 <= x && x < -0.375){
            value = -1.24277494216897627173 - 3.96044123930649416268 * x;
        } else if(0.375 <= x && x < 0.5 ){
            value = 1.11118460317702247406 - 2.51294751039500238824 * x;
        } else if(-1.0 <= x && x < -0.75){
            value = -0.491630451656180539739 - 0.491630451656180539739 * x;
        } else if(0.75 <= x && x < 1){
            value = 0.096859434680319146709 - 0.096859434680319146709 * x;
        } else if(0.5 <= x && x < 0.75) {
            value = -0.48429717340159573354 + 0.67801604276223402696 * x;
        } else if(-0.75 <= x && x < -0.5){
            value = 2.458152258280902698695 + 3.44141316159326377817 *x;
        } else if(-0.125 <= x && x < 0){
            value = 2.17124059336723766167 + 17.1233461124244526904 *x;
        } else if(-0.25 <= x && x < -0.125){
            value = 2.28083109759543704076 + 18.0000701462500477231 *x;
        } else if(0.125 <= x && x < 0.25){
            value = -4.78104753844255919661 + 22.6057320109684306418* x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0 <= x && x < 0.125) {
            value = -33.0125730435099442245;
        } else if(-0.375 <= x && x < -0.25) {
            value = -19.69261569230427143433;
        } else if(0.25 <= x && x < 0.375){
            value = -5.6124494201652150829;
        } else if(-0.5 <= x && x < -0.375){
            value = -3.96044123930649416268;
        } else if(0.375 <= x && x < 0.5 ){
            value = -2.51294751039500238824;
        } else if(-1.0 <= x && x < -0.75){
            value = -0.491630451656180539739;
        } else if(0.75 <= x && x < 1){
            value = -0.096859434680319146709;
        } else if(0.5 <= x && x < 0.75) {
            value =  0.67801604276223402696;
        } else if(-0.75 <= x && x <-0.5){
            value =  3.44141316159326377817;
        } else if(-0.125 <= x && x < 0){
            value = 17.1233461124244526904;
        } else if(-0.25 <= x && x < -0.125){
            value = 18.0000701462500477231;
        } else if(0.125 <= x && x < 0.25){
            value = 22.6057320109684306418;
        } else {
            value = 0.0;
        }
        
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_inner_evaluator1(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (0.125 <= x && x< 0.25) {
            value = 6.84500149587779311472 - 33.09125988795909139990* x;
        } else if(-0.375 <= x && x < -0.25){
            value = -4.757606755902492947821 - 13.01056641923786902976 *x;
        } else if(-0.5 <= x && x < -0.375){
            value = -1.761589978102878832866 - 5.02118834510556472321*x;
        } else if(0.5 <= x && x < 0.75){
            value = 0.953647148545425993185 - 1.335106007963596390459*x;
        } else if(-1.0 <= x && x < -0.75){
            value = -0.4993361296332690191592 - 0.4993361296332690191592*x;
        } else if(0.75 <= x && x < 1.0){
            value = -0.1907294297090851986370 + 0.1907294297090851986370*x;
        } else if(-0.125 <= x && x < 0.0){
            value = 5.797244214189106508582*x;
        } else if(0.0 <= x && x < 0.125){
            value = 21.66875207906325351789* x;
        } else if(-0.75 <= x && x < -0.5){
            value = 2.496680648166345095796 + 3.495352907432883134114*x;
        } else if(0.375 <= x && x < 0.5){
            value = -2.378803377951246473911 + 5.32979504502974854373*x;
        } else if(-0.25 <= x && x < -0.125){
            value = 0.0556540975457490632366 + 6.24247699455509901447*x;
        } else if(0.25 <= x && x < 0.375){
            value = -3.523179956205757665732 + 8.38146592037511172192*x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0.125 <= x && x< 0.25) {
            value = - 33.09125988795909139990;
        } else if(-0.375 <= x && x < -0.25){
            value =  - 13.01056641923786902976 ;
        } else if(-0.5 <= x && x < -0.375){
            value =  - 5.02118834510556472321;
        } else if(0.5 <= x && x < 0.75){
            value = - 1.335106007963596390459;
        } else if(-1.0 <= x && x < -0.75){
            value = - 0.4993361296332690191592;
        } else if(0.75 <= x && x < 1.0){
            value = 0.1907294297090851986370;
        } else if(-0.125 <= x && x <= 0.0){////ttttttttttest <=0.0
            value = 5.797244214189106508582;
        } else if(0.0 <= x && x < 0.125){
            value = 21.66875207906325351789;
        } else if(-0.75 <= x && x < -0.5){
            value =  3.495352907432883134114;
        } else if(0.375 <= x && x < 0.5){
            value =  5.32979504502974854373;
        } else if(-0.25 <= x && x < -0.125){
            value =  6.24247699455509901447;
        } else if(0.25 <= x && x < 0.375){
            value =  8.38146592037511172192;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_inner_evaluator2(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (0.625 <= x && x < 0.75) {
            value = 25.3256618915819297146 - 35.2183939040239842910*x;
        } else if(0.25 <= x && x < 0.375){
            value = 3.007148855493370156849 - 9.72312653910588252146*x;
        } else if(0.375 <= x && x < 0.5){
            value = 2.887238244145711420071 - 9.40336490884545922339*x;
        } else if(0.0 <= x && x < 0.125){
            value = -0.56347593702122728620*x;
        } else if (0.875 <= x && x < 1.0){
            value = -2.22797669417418531629 + 2.22797669417418531629 *x;
        } else if (0.125 <= x && x < 0.25){
            value = -0.717236204972206348034 + 5.17441370275642349806*x;
        } else if (0.75 <= x && x < 0.875){
            value = -5.94595223442177053803 + 6.47709159731428271256*x;
        } else if (0.5 <= x && x < 0.625){
            value = -22.32888385765284908917 + 41.0288792947516617951*x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0.625 <= x && x < 0.75) {
            value =  - 35.2183939040239842910;
        } else if(0.25 <= x && x < 0.375){
            value = - 9.72312653910588252146;
        } else if(0.375 <= x && x < 0.5){
            value = - 9.40336490884545922339;
        } else if(0.0 <= x && x < 0.125){
            value = -0.56347593702122728620;
        } else if (0.875 <= x && x < 1.0){
            value = 2.22797669417418531629 ;
        } else if (0.125 <= x && x < 0.25){
            value = 5.17441370275642349806;
        } else if (0.75 <= x && x < 0.875){
            value = 6.47709159731428271256;
        } else if (0.5 <= x && x < 0.625){
            value = 41.0288792947516617951;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_right_evaluator0(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (-0.375 <= x && x <-0.25) {
            value = -8.15586037317120966836 - 22.30372717533381322028*x;
        } else if (-0.5 <= x && x < -0.375){
            value = -3.01985486260715776817 - 8.60771248049634148643*x;
        } else if (-1.0 <= x && x < -0.75){
            value = -0.856000918427341983366 - 0.856000918427341983366*x;
        } else if (-0.125 <= x && x < 0){
            value = 9.93808794756765223095*x;
        } else if (-0.75 <= x && x < -0.5){
            value = 4.28000459213670991683 + 5.99200642899139388356*x;
        } else if (-0.25 <= x && x < -0.125){
            value = 0.0954065924458433055547 + 10.70134068713439867538*x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (-0.375 <= x && x <-0.25) {
            value = - 22.30372717533381322028;
        } else if (-0.5 <= x && x < -0.375){
            value =  - 8.60771248049634148643;
        } else if (-1.0 <= x && x < -0.75){
            value = - 0.856000918427341983366;
        } else if (-0.125 <= x && x < 0){
            value = 9.93808794756765223095;
        } else if (-0.75 <= x && x < -0.5){
            value =  5.99200642899139388356;
        } else if (-0.25 <= x && x < -0.125){
            value = 10.70134068713439867538;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_right_evaluator1(T x, unsigned short deriv)
{
    double value = 0.0;
    //x = power2(-1.0)*x;//anpassung
    if(deriv == 0){
        if(-0.375 <= x && x < -0.25){
            value = -4.2800045921367099168 - 11.984012857982787767*x;
        } else if(-0.5 <= x && x < -0.375){
            value = 0.8560009184273419834 + 1.712001836854683967*x;
        } else if (-0.125 <= x && x < 0.0){
            value = 3.7804491885886564460 + 19.876175895135304462*x;
        } else if (-0.25 <= x && x < -0.125){
            value = 3.8758557810344997515 + 20.639428634702050906*x;
        } else{
            value = 0.0;
        }
    } else if (deriv == 1) {
        if(-0.375 <= x && x < -0.25){
            value = - 11.984012857982787767;
        } else if(-0.5 <= x && x < -0.375){
            value =  1.712001836854683967;
        } else if (-0.125 <= x && x < 0.0){
            value =  19.876175895135304462;
        } else if (-0.25 <= x && x < -0.125){
            value =  20.639428634702050906;
        } else{
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

} // namespace lawa
