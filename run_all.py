import main as preproc
import compile_biokg
import build_benchmarks
import package_data

def main():
    print('Preprocessing data')
    preproc.preprocess_graph()
    print('Compiling graph')
    compile_biokg.compile_graph()
    print('Building benchmarks')
    build_benchmarks.main()
    print('Packaging data')
    package_data.package_all()

if __name__ == '__main__':
    main()