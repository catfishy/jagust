from collections import defaultdict
import pylab as P
import numpy as np
import itertools
from utils import *


def generate_Rousset_input(subjs, groupings, output_file):
    # write combinations
    outfile = open(output_file,'w')
    for a,b in itertools.product(subjs, groupings):
        subj_str = str(a).zfill(4)
        line = "%s;%s\n" % (subj_str,b)
        outfile.write(line)
    outfile.close()


def main(master_file):
    data = importMaster(master_file)
    # split by diagnosis at first av45 scan: 'Diag@AV45_long'
    by_diag = defaultdict(dict)
    for subj, d in data.iteritems():
        if int(subj) < 2000:
            continue
        diag = d['Diag@AV45_long'].strip()
        wcereb = d['AV45_wcereb'].strip()
        if diag == '':
            continue
        if wcereb == '':
            continue
        by_diag[diag][subj] = float(wcereb)

    AD = dict(by_diag['AD'])
    N = dict(by_diag['N'])

    # sort by whole cereb value: AV45_wcereb
    sorted_AD = sorted(AD.items(), key=lambda x: x[1])
    sorted_AD_vals = [_[1] for _ in sorted_AD]
    sorted_N = sorted(N.items(), key=lambda x: x[1])
    sorted_N_vals = [_[1] for _ in sorted_N]


    # print histograms
    P.figure()
    n, bins, patches = P.hist([sorted_AD_vals], 
                              70, 
                              histtype='bar',
                              label=['AD'])
    print "AD"
    for _ in zip(bins,n):
        print "\t%s" % (_,)
    P.legend()
    P.figure()
    n, bins, patches = P.hist([sorted_N_vals], 
                              70, 
                              histtype='bar',
                              label=['N'])
    print "N"
    for _ in zip(bins,n):
        print "\t%s" % (_,)

    P.legend()
    #P.show()

    # grab high, low, and around threshold (1.11)
    N_low = [(a,_) for a,_ in sorted_N if _ >= 1.0 and _ <= 1.03]
    N_mid = [(a,_) for a,_ in sorted_N if _ >= 1.09 and _ <= 1.13]
    N_high = [(a,_) for a,_ in sorted_N if _ >= 1.30 and _ <= 1.40]
    AD_low = [(a,_) for a,_ in sorted_AD if _ >= 0.9 and _ <= 1.03]
    AD_mid = [(a,_) for a,_ in sorted_AD if _ >= 1.09 and _ <= 1.13]
    AD_high = [(a,_) for a,_ in sorted_AD if _ >= 1.35 and _ <= 1.42]

    print "N_LOW"
    for _ in N_low: print _
    print "N_MID"
    for _ in N_mid: print _
    print "N_HIGH"
    for _ in N_high: print _
    print "AD_LOW"
    for _ in AD_low: print _
    print "AD_MID"
    for _ in AD_mid: print _
    print "AD_HIGH"
    for _ in AD_high: print _

    # return list of all subjects
    all_AD = [_[0] for _ in itertools.chain(AD_low, AD_mid, AD_high)]
    all_N = [_[0] for _ in itertools.chain(N_low, N_mid, N_high)]

    print "%s AD and %s N" % (len(all_AD), len(all_N))
    return (all_N, all_AD)

if __name__ == "__main__":
    '''
    # To find low/mid/high scans to run PVC on
    master = "../FDG_AV45_COGdata_07_28_15.csv"
    all_N, all_AD = main(master)
    print "AD: %s" % all_AD
    print "N: %s" % all_N

    subj = all_N + all_AD
    print subj
    sys.exit(1)
    '''
    # To generate the Rousset SGE input file
    subj_list = ['003-S-0907', '003-S-4152', '011-S-4827', '003-S-2374', '011-S-0021', '011-S-4845', '003-S-4136', '003-S-4081', '011-S-4893', '003-S-1057', '011-S-4906', '003-S-4119', '011-S-4912', '003-S-4350', 
                 '003-S-0908', '003-S-4288', '011-S-4949', '003-S-0981', '003-S-1122', '003-S-1074', '002-S-0413', '009-S-4337', '009-S-4612', '009-S-5027', '011-S-0023', '012-S-4012', '002-S-0685', '009-S-5037', '011-S-2274', '012-S-4026', '002-S-1155', '009-S-4741', '009-S-5125', '011-S-4075', '012-S-4094', 
                 '002-S-0729', '011-S-4105', '012-S-4128', '003-S-4555', '003-S-4354', '003-S-4441', '012-S-4188', '003-S-4644', '005-S-0546', '009-S-5147', '011-S-4120', '012-S-4545', '005-S-4168', '005-S-2390', 
                 '005-S-0610', '005-S-4185', '005-S-4707', '006-S-1130', '006-S-0498', '002-S-1261', '005-S-0602', '009-S-4359', '009-S-4814', '009-S-5176', '011-S-4222', '012-S-4643', '006-S-4546', '007-S-4272', '007-S-4488', '007-S-2058', '007-S-2394', '007-S-2106', '007-S-4387', '007-S-4516', '007-S-4467', 
                 '007-S-1206', '007-S-4568', '009-S-5224', '011-S-4235', '012-S-4849', '007-S-4611', '012-S-4987', '007-S-0101', '012-S-5121', '007-S-4620', '007-S-4637', '012-S-5157', '009-S-0842', '009-S-1030', '009-S-2381', '002-S-4654', '002-S-1268', '005-S-0553', '009-S-4388', '009-S-4903', '009-S-5252', 
                 '011-S-4278', '012-S-5195', '010-S-0419', '941-S-4377', '003-S-4373', '002-S-5256', '003-S-5187', '006-S-4449', '006-S-4679', '006-S-4515', '006-S-4150', '006-S-4153', '006-S-4346', '006-S-4357', '002-S-5230', '006-S-4192', '002-S-4213', '002-S-2043', '002-S-4219', '002-S-4171', '002-S-2010', 
                 '002-S-5178', '002-S-4251', '002-S-2073', '002-S-4270', '002-S-5018', '002-S-4521', '002-S-4746', '002-S-4799', '002-S-4225', '002-S-4237', '002-S-4447', '009-S-4324', '007-S-0698', '006-S-4867', '006-S-4960', '006-S-5153', '006-S-4713', '009-S-2208', '009-S-0751', '002-S-4229', '002-S-1280', 
                 '002-S-4262', '006-S-4485', '009-S-4530', '009-S-4958', '010-S-0420', '011-S-4366', '012-S-5213', '014-S-4668', '016-S-0702', '016-S-1117', '016-S-1326', '016-S-2007', '016-S-2031', '016-S-2284', '016-S-4009', '016-S-4097', '016-S-4121', '016-S-4353', '016-S-4583', '016-S-4584', '016-S-4601', 
                 '016-S-4591', '016-S-4638', '016-S-4646', '016-S-4688', '002-S-0295', '002-S-4473', '006-S-4363', '009-S-4543', '009-S-5000', '010-S-4345', '011-S-4547', '013-S-5171', '014-S-4039', '014-S-4328', '018-S-2133', '018-S-4400', '019-S-4252', 
                 '019-S-4835', '021-S-4924', '022-S-4291', '023-S-2068', '023-S-4448', '021-S-0159', '023-S-4501', '021-S-0984', '021-S-2077', '021-S-2100', '021-S-2124', '021-S-2125', '021-S-2142', '021-S-2150', '021-S-4245', '021-S-4254', '021-S-4276', '021-S-4335', '021-S-4402', '021-S-4419', '021-S-4421', 
                 '021-S-4558', '021-S-4659', '021-S-4718', '021-S-4744', '021-S-4857', '022-S-0096', '022-S-0130', '014-S-4058', '014-S-4401', '018-S-2138', '018-S-4597', '019-S-4285', '019-S-5012', '022-S-2087', '022-S-4320', '023-S-4020', '023-S-4502', '023-S-0042', '023-S-0031', '023-S-0058', '023-S-0061', 
                 '023-S-0126', '023-S-0331', '023-S-4034', '023-S-4796', '023-S-0926', '023-S-1046', '023-S-1190', '013-S-4985', '014-S-0658', '014-S-4079', '014-S-4576', '018-S-2155', '018-S-4696', '019-S-4293', '019-S-5019', '022-S-2167', '022-S-4444', '023-S-4035', '023-S-5120', '024-S-0985', '024-S-1063', '018-S-2180', 
                 '018-S-4809', '019-S-4367', '019-S-5242', '022-S-2263', '022-S-4805', '023-S-4115', '023-S-5241', '027-S-0074', '027-S-0118', '027-S-0120', '027-S-0307', '027-S-0408', '027-S-0644', '027-S-0835', '027-S-1045', '027-S-2183', '027-S-2219', '027-S-2245', '027-S-2336', '027-S-4729', '027-S-4757', 
                 '027-S-4801', '027-S-4802', '027-S-4804', '029-S-0845', '029-S-1218', '029-S-1318', '029-S-2376', '029-S-2395', '029-S-4279', '029-S-4290', '029-S-4307', '029-S-4327', '029-S-4384', '029-S-4385', '029-S-4585', '029-S-4652', '031-S-0618', '031-S-0830', '031-S-0867', '013-S-1186', '013-S-2324', 
                 '013-S-2389', '013-S-4395', '013-S-4579', '013-S-4616', '013-S-4917', '013-S-5071', '014-S-2185', '014-S-4080', '014-S-4577', '018-S-4313', '018-S-4868', '019-S-4477', '020-S-4920', '022-S-2379', '022-S-4922', '023-S-4122', '024-S-2239', '032-S-0214', '032-S-0479', '014-S-4093', '014-S-4615', '018-S-4349', 
                 '018-S-4889', '019-S-4548', '020-S-5140', '022-S-4173', '022-S-5004', '023-S-4164', '024-S-4084', '033-S-0741', '033-S-0906', '033-S-0920', '033-S-0923', '033-S-1016', '019-S-4549', '020-S-5203', '022-S-4196', '023-S-0376', '023-S-4241', '024-S-4158', '035-S-0292', '035-S-0997', '013-S-4268', '013-S-4580', '013-S-4595', '013-S-4791', '013-S-5137', '014-S-2308', '014-S-4263', '016-S-5251', '018-S-4399', '018-S-5250', '019-S-4680', '021-S-0626', '022-S-4266', '023-S-0887', '023-S-4243', 
                 '024-S-4169', '036-S-4820', '036-S-5271', '037-S-4302', '037-S-5162', '041-S-4200', '041-S-4989', '037-S-0150', '037-S-0454', '037-S-0552', '033-S-4176', '033-S-5013', '035-S-2061', '035-S-4256', '035-S-4785', '036-S-4538', '036-S-4878', '036-S-5283', '037-S-4308', '037-S-5222', '041-S-4271', '041-S-5026', '041-S-0679', '041-S-1010', 
                 '041-S-1418', '041-S-1425', '032-S-2240', '032-S-4386', '032-S-4921', '033-S-4177', '033-S-5017', '035-S-2074', '035-S-4414', '036-S-2378', '036-S-4562', '036-S-4894', '037-S-4001', '037-S-4381', '041-S-4004', '041-S-4427', '041-S-5078', '052-S-0671', '052-S-0989', '052-S-1346', '052-S-1352', '052-S-2249', 
                 '052-S-4626', '053-S-0919', '037-S-4410', '041-S-4014', '041-S-4510', '041-S-5082', '057-S-0934', '057-S-1007', '057-S-1269', '041-S-5097', '057-S-4897', '067-S-0056', '067-S-0059', '033-S-5087', '035-S-2199', '035-S-4464', '036-S-2380', '036-S-4714', '036-S-4899', '037-S-4015', '037-S-4706', '041-S-4037', '041-S-4513', '041-S-5100', '068-S-0127', '068-S-0872', '036-S-5063', '037-S-4028', '037-S-4750', '041-S-4041', '041-S-4629', '041-S-5131', '068-S-2315', '031-S-4947', '032-S-2247', 
                 '032-S-4429', '032-S-5263', '033-S-4179', '033-S-5198', '035-S-4082', '035-S-4582', '036-S-4389', '036-S-4715', '036-S-5112', '037-S-4030', '037-S-4770', '041-S-4051', '041-S-4720', '041-S-5141', '072-S-0315', '041-S-5204', '072-S-1380', '041-S-4874', '041-S-5244', '072-S-2037', '041-S-4060', '041-S-4876', '041-S-5253', '072-S-2093', '024-S-4223', '024-S-4280', '024-S-4392', '024-S-4674', '024-S-4905', 
                 '024-S-5054', '024-S-5290', '031-S-1066', '031-S-2018', '031-S-2022', '031-S-2233', '031-S-4005', '031-S-4021', '031-S-4024', '031-S-4029', '031-S-4032', '031-S-4149', '031-S-4203', '031-S-4474', '031-S-4496', '031-S-4721', '032-S-1169', '032-S-4277', '032-S-4755', '032-S-5289', '033-S-4505', '033-S-5235', '035-S-4085', '035-S-4783', '036-S-4430', '036-S-4736', '036-S-5210', '037-S-4146', '037-S-4879', '041-S-4138', '041-S-4877', '051-S-4929', '073-S-4312', '031-S-4042', '031-S-4194', 
                 '031-S-4218', '031-S-4476', '031-S-4590', '032-S-2119', '032-S-4348', '032-S-4823', '033-S-1098', '033-S-4508', '033-S-5259', '035-S-4114', '035-S-4784', '036-S-4491', '036-S-4740', '036-S-5248', '037-S-4214', '037-S-5126', '041-S-4143', '041-S-4974', '051-S-4980', '094-S-2201', '094-S-2216', '094-S-2238', '094-S-2367', '094-S-4089', 
                 '094-S-4162', '094-S-4234', '094-S-4282', '094-S-4434', '094-S-4503', '094-S-4560', '072-S-2164', '094-S-4649', '094-S-4737', '098-S-0160', '098-S-0171', '098-S-0172', '098-S-0667', '098-S-0896', '098-S-2047', '098-S-2052', '098-S-2079', '098-S-4003', '098-S-4018', '098-S-4050', '098-S-4059', 
                 '098-S-4201', '098-S-4215', '098-S-4275', '098-S-4506', '099-S-0051', '099-S-0291', '099-S-0352', '053-S-5272', '053-S-5296', '067-S-2195', '067-S-2301', '067-S-4054', '067-S-4184', '067-S-4310', '067-S-4767', '067-S-5159', '067-S-5212', '068-S-2184', '068-S-2248', '068-S-4217', '068-S-4859', '070-S-4856', '072-S-4007', '100-S-0047', '100-S-0069', '100-S-0296', '100-S-1226', '100-S-1286', '068-S-4968', '070-S-5040', '072-S-4057', '109-S-2200', '109-S-2237', '109-S-2278', '109-S-4380', 
                 '109-S-4455', '109-S-4499', '109-S-4531', '109-S-4594', '114-S-0173', '114-S-0378', '114-S-0416', '114-S-1106', '114-S-1118', '072-S-2026', '072-S-4063', '116-S-0361', '116-S-0657', '116-S-1232', '116-S-4010', '116-S-4043', '068-S-2316', '068-S-4274', '068-S-5146', '072-S-2027', '072-S-4102', '116-S-4209', '068-S-2187', '068-S-4061', 
                 '068-S-4332', '070-S-4692', '072-S-2070', '072-S-4103', '123-S-0106', '123-S-0113', '123-S-1300', '067-S-5160', '068-S-2168', '068-S-2193', '068-S-4067', '068-S-4340', '070-S-4708', '072-S-2072', '072-S-4131', '126-S-0680', '126-S-0709', '126-S-1187', '126-S-2360', '126-S-2405', '126-S-2407', '126-S-4458', '126-S-4494', '126-S-4507', '126-S-4514', '126-S-4675', '126-S-4686', '126-S-4712', '127-S-0112', '127-S-0260', '127-S-0925', '127-S-1419', '127-S-1427', '127-S-2213', '127-S-2234', 
                 '127-S-4148', '127-S-4197', '127-S-4198', '127-S-4210', '127-S-4240', '127-S-4301', '127-S-4500', '127-S-4604', '127-S-4624', '127-S-4645', '127-S-4765', '127-S-4843', '128-S-0135', '128-S-0200', '128-S-0225', '128-S-0227', '128-S-0229', '128-S-0230', '128-S-0272', '128-S-0545', '128-S-0863', '128-S-1043', '128-S-1408', '051-S-5005', 
                 '051-S-5285', '053-S-2357', '053-S-2396', '053-S-4557', '053-S-4578', '053-S-4661', '053-S-4813', '053-S-5070', '053-S-5202', '053-S-5208', '053-S-5287', '057-S-2398', '067-S-2196', '067-S-2304', '067-S-4072', '067-S-4212', '067-S-4728', '067-S-4782', '067-S-5205', '068-S-2171', '068-S-2194', '068-S-4134', '068-S-4424', '070-S-4719', '072-S-2083', '072-S-4206', '129-S-1246', '129-S-2332', '068-S-4174', '068-S-4431', '070-S-4793', '072-S-2116', '072-S-4226', '100-S-5106', '116-S-4453', 
                 '123-S-4780', '130-S-0232', '130-S-0285', '130-S-0289', '130-S-0886', '073-S-2182', '073-S-2225', '073-S-4216', '073-S-4311', '073-S-4393', '073-S-4540', '073-S-4739', '073-S-4825', '073-S-5227', '082-S-4090', '082-S-4339', '082-S-5184', '099-S-2042', '099-S-4076', '099-S-4480', '100-S-5280', '116-S-4483', '123-S-4806', '131-S-0123', '131-S-0441', '073-S-4552', '073-S-4762', '073-S-4853', '082-S-2099', '082-S-4208', '082-S-4428', '082-S-5279', '099-S-2063', '099-S-4086', '099-S-4498', 
                 '114-S-4379', '116-S-4625', '123-S-4904', '136-S-0107', '136-S-0186', '136-S-0873', '099-S-4104', '099-S-4565', '114-S-4404', '116-S-4635', '127-S-5228', '137-S-0301', '137-S-0722', '137-S-0800', '137-S-0972', '137-S-0973', '137-S-0994', '137-S-1414', '073-S-2190', '073-S-2264', '073-S-4259', '073-S-4360', '073-S-4403', '073-S-4559', '073-S-4777', '073-S-5023', '082-S-2121', '082-S-4224', '082-S-5014', '082-S-5282', '099-S-2146', '099-S-4157', '099-S-4994', '114-S-5047', '116-S-4732', 
                 '127-S-5266', '141-S-0717', '128-S-2002', '141-S-1004', '141-S-1052', '141-S-1255', '141-S-1378', '072-S-4383', '072-S-4390', '072-S-4391', '072-S-4394', '072-S-4445', '072-S-4462', 
                 '072-S-4465', '072-S-4522', '072-S-4539', '072-S-4610', '072-S-4613', '072-S-4694', '072-S-4769', '072-S-4871', '072-S-4941', '072-S-5207', '073-S-0089', '073-S-0311', '073-S-0746', '073-S-2153', '073-S-2191', '073-S-4155', '073-S-4300', '073-S-4382', '073-S-4443', '073-S-4614', '073-S-4795', '073-S-5167', '082-S-2307', '082-S-4244', '082-S-5029', '094-S-4630', '099-S-2205', '099-S-4202', '100-S-4469', '114-S-5234', '116-S-4855', '128-S-2003', '018-S-0055', '018-S-0142', '099-S-4205', 
                 '100-S-4512', '116-S-4092', '116-S-4898', '128-S-2011', '052-S-4885', '052-S-4807', '123-S-2055', '128-S-2036', '057-S-4888', '100-S-4556', '116-S-4167', '123-S-2363', '128-S-2045', '114-S-0166', '128-S-2057', '126-S-4743', '126-S-4891', '126-S-4896', '127-S-0259', '127-S-4844', '100-S-5075', '116-S-4175', '123-S-4096', '128-S-2123', '021-S-0276', '130-S-0505', '003-S-4872', '003-S-4892', '003-S-4900', '007-S-4911', '099-S-4463', '100-S-5091', '116-S-4195', '123-S-4127', '128-S-2130', 
                 '016-S-4887', '016-S-4902', '099-S-4022', '099-S-4475', '100-S-5096', '116-S-4199', '123-S-4170', '128-S-2151', '027-S-4869', '027-S-4873', '027-S-4919', '100-S-5102', '116-S-4338', '123-S-4526', '128-S-2220', '130-S-4982', '135-S-4676', '137-S-4466', '057-S-4909', '130-S-4352', '130-S-4984', '135-S-4689', '137-S-4482', '127-S-4928', '128-S-4599', '128-S-4636', '128-S-4774', '129-S-4287', '130-S-4405', '130-S-4990', '135-S-4722', '137-S-4520', 
                 '067-S-4918', '127-S-4940', '135-S-4723', '137-S-4536', '005-S-4910', '005-S-5038', '128-S-4792', '129-S-4369', '130-S-4415', '130-S-4997', '135-S-4863', '137-S-4587', '016-S-4951', '016-S-4952', '016-S-4963', 
                 '016-S-5007', '128-S-4832', '129-S-4371', '130-S-4417', '130-S-5006', '135-S-4954', '137-S-4596', '027-S-4926', '027-S-4936', '027-S-4938', '027-S-4943', '027-S-4955', '027-S-4962', '027-S-4966', '130-S-5059', '135-S-5015', '137-S-4623', '027-S-4964', '128-S-4653', '128-S-4842', '129-S-4396', '130-S-4468', '130-S-5142', '135-S-5113', '137-S-4631', '073-S-4986', '130-S-5175', '135-S-5269', '137-S-4632', '094-S-4858', '135-S-5273', '137-S-4672', '127-S-4992', '127-S-5028', '127-S-5056', 
                 '128-S-4607', '128-S-4671', '128-S-5066', '129-S-4422', '130-S-4542', '130-S-5231', '135-S-5275', '137-S-4678', '052-S-4944', '052-S-4945', '130-S-5258', '136-S-4189', '137-S-4756', 
                 '127-S-5058', '127-S-5067', '128-S-2314', '128-S-4553', '128-S-4571', '128-S-4586', '128-S-4609', '128-S-4742', '129-S-0778', '130-S-2373', '130-S-4589', '131-S-5138', '136-S-4269', '137-S-4815', '003-S-5130', '005-S-5119', '130-S-2391', '130-S-4605', '135-S-4281', '136-S-4408', '137-S-4816', '016-S-5031', '016-S-5057', '021-S-5099', '021-S-5129', '027-S-5093', '027-S-5109', '027-S-5110', '027-S-5118', '027-S-5127', '029-S-5135', '128-S-4745', '129-S-2347', '130-S-2403', '130-S-4641', 
                 '135-S-4309', '136-S-4433', '137-S-4852', '052-S-4959', '052-S-5062', '135-S-4356', '136-S-4517', 
                 '137-S-4862', '127-S-5095', '127-S-5132', '127-S-5185', '130-S-4660', '135-S-4406', '137-S-4211', '141-S-0767', '007-S-5196', '130-S-4730', '135-S-4446', '137-S-4258', '141-S-2210', '021-S-5177', '021-S-5194', '129-S-4073', '130-S-4250', '130-S-4817', '135-S-4489', '137-S-4299', '141-S-2333', '057-S-5199', '128-S-4772', '129-S-4220', '130-S-4294', '130-S-4883', '135-S-4566', '137-S-4303', '141-S-4053', '126-S-5214', '141-S-4160', '029-S-5166', '141-S-4232', '003-S-5154', '130-S-4343', 
                 '130-S-4925', '135-S-4598', '137-S-4331', '141-S-4426', '003-S-5165', '141-S-4438', '003-S-5209', '007-S-5265', '130-S-4971', '135-S-4657', '137-S-4351', '141-S-4456', '021-S-5236', '021-S-5237', '027-S-5169', '027-S-5197', '029-S-5158', '029-S-5219', 
                 '127-S-5200', '127-S-5218', '018-S-5240', '027-S-5277', '126-S-5243', '027-S-5170', '027-S-5288', '057-S-5292', '141-S-4711', '141-S-4803', '141-S-4907', '141-S-4976', '153-S-2109', '153-S-2148', '153-S-4077', '153-S-4125', '153-S-4133', '153-S-4139', '153-S-4151', '153-S-4159', '153-S-4172', '153-S-4297', '153-S-4372', '153-S-4621', '153-S-4838', '153-S-5261', '153-S-5267', '941-S-1195', '941-S-1202', '941-S-2060', '941-S-4036', '941-S-4066', '941-S-4100', '941-S-4187', '941-S-4255', 
                 '941-S-4292', '941-S-4365', '941-S-4376', '941-S-4420', '941-S-4764', '941-S-5124', '941-S-5193', '027-S-5083', '027-S-5079', '114-S-2392', '137-S-0668', '128-S-5123', '021-S-0337']
    subj = [int(_.split('-')[-1]) for _ in subj_list]
    groupings = ['/home/jagust/ahorng/matlab/pvc/groupings_Scan2/grouping_2.mat',
                 '/home/jagust/ahorng/matlab/pvc/groupings_Scan2/grouping_4.mat',
                 '/home/jagust/ahorng/matlab/pvc/groupings_Scan2/grouping_agghigh.mat',
                 '/home/jagust/ahorng/matlab/pvc/groupings_Scan2/grouping_agglowtwo.mat']
    output = 'input_scan2.txt'
    generate_Rousset_input(subj, groupings, output)

    groupings = ['/home/jagust/ahorng/matlab/pvc/groupings_BL/grouping_2.mat',
                 '/home/jagust/ahorng/matlab/pvc/groupings_BL/grouping_4.mat',
                 '/home/jagust/ahorng/matlab/pvc/groupings_BL/grouping_agghigh.mat',
                 '/home/jagust/ahorng/matlab/pvc/groupings_BL/grouping_agglowtwo.mat']
    output = 'input_bl.txt'
    generate_Rousset_input(subj, groupings, output)

    groupings = ['/home/jagust/ahorng/matlab/pvc/groupings_Scan3/grouping_2.mat',
                 '/home/jagust/ahorng/matlab/pvc/groupings_Scan3/grouping_4.mat',
                 '/home/jagust/ahorng/matlab/pvc/groupings_Scan3/grouping_agghigh.mat',
                 '/home/jagust/ahorng/matlab/pvc/groupings_Scan3/grouping_agglowtwo.mat']
    output = 'input_scan3.txt'
    
    generate_Rousset_input(subj, groupings, output)
    sys.exit(1)

    '''
    output = 'input_agglow.txt'
    generate_Rousset_input(subj, ['/home/jagust/ahorng/matlab/pvc/groupings_082715/grouping_agglow.mat'], output)
    '''

    '''
    error_subj = [4022, 4376, 4309, 1202, 5177, 4730, 5261, 4963, 337, 4657, 4621, 4553, 4269, 4356, 4351, 4815, 4420, 4426, 4100, 5132, 5130, 4187, 4907, 5123, 4172, 5119, 4772, 4589, 4077, 4073, 5165, 4294, 5057, 1195, 4258, 4817, 4250, 4255, 4918, 2060, 5170, 5079, 4220, 4292, 4408, 5038, 5267, 2148, 4036, 5031, 4139, 4133, 4433, 4722, 4542, 668, 5236, 5237, 4343, 4438, 4463, 5169, 5109, 55, 4862, 5292, 4660, 2151, 2392, 5214, 5219, 5199, 4160, 5197, 2333, 5093, 4303, 4959, 4151, 4571, 4125, 5067, 5062, 5243, 4232, 4406, 4609, 4976, 778, 4971, 5218, 5099, 4127, 5095, 5288, 5265, 2347, 4925, 4053, 4711, 4159, 4372, 2205, 5196, 5185, 5110, 4446, 4371, 4883, 4816, 4838, 2109, 4764, 5118, 5240, 4598, 2403, 5200, 5209, 4331, 5193, 5194, 4517, 5277, 4211, 4566, 4641, 2314, 5138, 5135, 4803, 4417, 4900, 2391, 5129, 5124, 5127, 4299, 5083, 4297, 4489, 5158, 4066, 4742, 4745, 5154, 4605, 2210, 767, 2373, 4365, 4456, 2130, 4990, 5166, 4281, 4852, 4586]
    output = 'input_errors.txt'
    generate_Rousset_input(error_subj, groupings, output)
    '''