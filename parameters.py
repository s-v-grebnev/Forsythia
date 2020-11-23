"""
(c) 2020 Sergey Grebnev, s.v.grebnev@yandex.ru
"""
params = {
    #
    # Parameters from C. Costello's eprint 2019/1321
    #
    'costello': {'lA': 2, 'lB': 3, 'eA': 4, 'eB': 3, 'f': 1, 'A': [423, 329], 'C': [1, 0],
                 'xp2': [248, 100], 'yp2': [199, 304], 'zp2': [1, 0],
                 'xq2': [394, 426], 'yq2': [79, 51], 'zq2': [1, 0],
                 'xp3': [275, 358], 'yp3': [104, 410], 'zp3': [1, 0],
                 'xq3': [185, 20], 'yq3': [239, 281], 'zq3': [1, 0]
                 },

    #
    # Same as above, but starting curve is D. Koshelev's E_19
    #
    'toy_example': {'lA': 2, 'lB': 3, 'eA': 4, 'eB': 3, 'f': 1, 'A': [104, 0], 'C': [1, 0],
                    'xp2': [361, 304], 'yp2': [157, 414], 'zp2': [299, 113],
                    'xq2': [34, 368], 'yq2': [103, 124], 'zq2': [203, 212],
                    'xp3': [176, 172], 'yp3': [106, 319], 'zp3': [0, 70],
                    'xq3': [26, 154], 'yq3': [281, 18], 'zq3': [416, 298]
                    },
    #
    # Large-scale test
    #
    'p485': {'lA': 2, 'lB': 3, 'eA': 242, 'eB': 152, 'f': 5,
             'A': [
                 12496973461376851434681477318484969222675795134827832485840725601663384787637413718208157396895460518243784563742278285794546627758175340527396825,
                 0],
             'C': [1, 0],

             'xp2': [
                 61376322485193188540753756261850807725490692389731134988448397485321223953551161389388954517165523301620173196993441012542990173519935671131827257,
                 113151541168782510447152323635947283584254776788058488834564507058893980642337983486838759880014906863247921335846729412052486481124279441065058713],
             'yp2': [
                 107663114781372711691553296029737605188059195482326583207348227009500608879620657995495746524346759760624570984943202402251912952274383696557773302,
                 95791431627200553537776218364588635907583658118955468156462283620017536966324426045366054337441034018206938679120236614406204941034734783124741407
             ],
             'zp2': [
                 1058610704548972272547381304767672486462742330037705550024400602674372870824713998035233735808010856966768029869424703160796868878813369568350349,
                 95663897601708622428369826147128976160466116783115917587212759758956355646614374672705738162397522681492094197609402773078696304735040671588228405
             ],
             'xq2': [
                 64226903097778179984070510878577589775431536199692577851994197484986738666006152675639779142565915696337062761064811457293530761615107948099236735,
                 38110065820820119259072438212942338594066660820454172385667222365385088192058161071849371790209396277761158948013089762074214395625843991012901224
             ],
             'yq2': [
                 76581930107679640466716332083396778994046245863658537337701291399701582484518225178892892905418608287615132197560860002784009979931352046822930363,
                 21464785142317824123083304988389872008099660024808394917433432992755991326232710096058346519318908045230095879832002119070020857530599524484702639],
             'zq2': [
                 55254424798345243979201491759690485653998388700380352496419432797144733811820537063963294783185724778530534524522357233929868792733936087959297836,
                 93526513507719082431634542918788862313879735403740049239596795048923776478568013439683094270752002123543651977522387411639678931823556397541379744],
             'xp3': [
                 88368389078295434508330208033010337761324388277765111735918172956538284321312068238245323009703500562212937796246327757713573351393057936783484023,
                 105744384942355688522859515375450839555115012935291023378854985859031179379887206587854338983848883792571907274225505420984555928802320741554311780],
             'yp3': [
                 38675174826835071770590073799706381074249907240232339162723201682586765458813636046633961449681997524567717612123922821446991255517053044234770111,
                 10022851622984190757257737960862238711366173580717733682765669801648912451203940009183928618559807186883590896717174185528262159571156702933929722],
             'zp3': [
                 67864368399181361155125252409426024530428739067947021791253287941492953287847997550195240594216117178250383842686463597590562173719885766701168032,
                 55796620791380653242703974561229687526997479297786108306696542058884586338073983557707499235076402750922453975830621197143831176894170887315770910],
             'xq3': [
                 83136069274649107397411910270099411864305577923718255289759464671223574535823374925254900363542325208625593188743312791408937711829234296149257376,
                 103222198997057791264217549626741420352854212672901349031557582268686319615058029357930595113048679430659931480683289789015263383711245401297138098],
             'yq3': [
                 1490584805919108709438838603222656387518537507388450514632828203412939243741169214278051127021401171627074827135085454366174067053163737598521694,
                 51659041050386386530527849766968994972278275811368912362947342218886793073132134592145653027963586586503132033115594959302574496119137001879436143],
             'zq3': [
                 15661197771654125795760151382296793949619438973156538619717525192114724572983717763380698540351996603899573106609614881627372815305692894537055868,
                 43570449862743614870748251652845615611604237525954885195348158475641366950809972885319911843753451832357616167854406765973402954521417177382819764]
             },
    'forsythia80': {'lA': 2, 'lB': 3, 'eA': 137, 'eB': 84, 'f': 1,
                    'A': [1697515660496260333664277950769004218366836292978277269037599020171750501119102405,
                        0],
                    'C': [1, 0],
    'xp2': [275054030742979557334472621739602124103041806982787060420183725308440589518956919,
            1197135460420307614825932378602332597171569965751622984352657305089633657273652277],
    'yp2': [850241066313357956000226093312616073708173909912271684766421784986680021655693717,
            1495027069330966568472240789545176489057626594566039984682225976374807036776099748],
    'zp2': [1558443940022337609329074033043754990259683811300904869155295865231693971532567047,
            415030868081208009808755776461888553346801439135405426105642235753028632801132971],
    'xq2': [684209869270106587359151829681501020811363343390393290847601448142956880576003004,
            1487392302178946120337364613883142729563168332804489918060135924880019431445391942],
    'yq2': [404029536748315088250722663846834101316577173530257439370670546449432833652863969,
            1790160113049455324983290956786119033179004146300393653050373458448097788194785376],
    'zq2': [725426772614053576913875824809966048521506075201977583213216482040798098766053700,
            133162539602818633576403006442460143165927428993380624413067471587054584530661578],
    'xp3': [1874501549628802227094010464821792720079597264714503778336531897093381721318231101,
            123547905020588398058816483168196821015330776577490158717505544277834994271300495],
    'yp3': [187838296951638101708179236469019302111621334628939352985257112022743383774927514,
            1512542006319805574752928797919803324395390189691336891138606141267521949216470580],
    'zp3': [1779022630404946457965055242334369693167398755703772890237143979646138270753257174,
            1531026453823957401679644451849813350363453262037338368957042708023475326777130385],
    'xq3': [1459218276066324788844669144152061914950504413752532397867055750258839873039135526,
            1641888668580189721444640150434291222084957036284832174036148053878523395844063189],
    'yq3': [993842882207021781600129337812749372273039522369063106762379319095037584771278822,
            1673483427444379450811966194148240047262342026061987220201119910660141782912831561],
    'zq3': [1756247221970242619933516828558211142817228548132681730438431614826301215892652673,
            2033180697801576059018957017489135843212785691829562911577277220173335725225571480],
    }

}

