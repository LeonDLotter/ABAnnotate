
def get_aba_data(atlas, save_dir, save_prefix, **kwargs):

    from os.path import join
    try:
        from abagen import get_expression_data
    except:
        print("Error importing abagen toolbox. To install, run 'pip install abagen' in terminal.")

    # get roi-wise mRNA expression data (if not already done so, will 
    # download complete human ABA data @ ~4GB) through the abagen toolbox
    print('Loading ABA data. This will take some time. If not done before, downloading complete human ABA dataset (~4 GB).')
    aba_parcellated, method = get_expression_data(atlas=atlas, return_report=True, **kwargs)
    print(f'{aba_parcellated.shape[1]} genes in processed AHBA dataset.')

    # save
    aba_parcellated_file = join(save_dir, save_prefix+'_data.csv')
    report_file = join(save_dir, save_prefix+'_report.md')

    print(f'Saving parcel-wise data to: {aba_parcellated_file}')
    aba_parcellated.to_csv(aba_parcellated_file)

    print(f'Saving method report to: {report_file}')
    t = open(report_file, 'w')
    t.write(method)
    t.close()

    print('Finished importing ABA data.')

