import argparse
import sys

import interactions

def main():
    parser = argparse.ArgumentParser(description='CLI for calculating skyglow')
    parser.add_argument('action', type=str, help='sgmap_single/kernel_lib/sgmap_multiple')
    parser.add_argument('-t', '--target', dest='target', type=int, help='Gene to run query on (entrez id)')
    parser.add_argument('-o', '--order', dest='order', type=int, default=1, help='Order of interacting genes away from target')

    args = parser.parse_args()
    if args.operation != 'populate' and args.operation != 'query':
        parser.print_usage(sys.stderr)
        print 'Error: operation must be populate or query'
    if args.operation == 'query' and args.target == None:
        parser.print_usage(sys.stderr)
        print 'Error: target must be specified for query operation'

    if args.operation == 'populate':
        interactions.populate()
    elif args.operation == 'query':
        result = interactions.query(args.target, args.order)
        msg = '{}-order interacting genes with gene {}:'.format(args.order, args.target)
        genes = []
        for d in result:
            genes.append(d['n_order_gene.id'])
        print msg, ', '.join(map(str,genes))


if __name__ == '__main__':
    main()
