import os
from Bio import SeqIO
from typing import List, Dict


class Node:
    def __init__(self, kmer: str):
        self.kmer = kmer
        self.in_edges = []
        self.out_edges = []


class Edge:
    def __init__(self, sequence: str):
        self.sequence = sequence  # k+1 мер (или длиннее после сжатия)
        self.coverage = 1.0  # Глубина покрытия (float, так как будет усредняться)
        self.source: Node = None
        self.target: Node = None


class DeBruijnGraph:
    def __init__(self, k: int):
        self.k = k
        self.nodes: Dict[str, Node] = {}
        self.edges: List[Edge] = []

    def get_or_create_node(self, kmer: str) -> Node:
        if kmer not in self.nodes:
            self.nodes[kmer] = Node(kmer)
        return self.nodes[kmer]

    def add_sequence(self, sequence: str):
        """Разбивает последовательность на k-меры и строит (k+1)-ребра."""
        # Проходим скользящим окном размером k+1
        for i in range(len(sequence) - self.k):
            left_kmer = sequence[i: i + self.k]
            right_kmer = sequence[i + 1: i + self.k + 1]
            edge_seq = sequence[i: i + self.k + 1]

            node_left = self.get_or_create_node(left_kmer)
            node_right = self.get_or_create_node(right_kmer)

            # Проверяем, существует ли уже такое ребро между этими узлами
            existing_edge = None
            for e in node_left.out_edges:
                if e.target == node_right:
                    existing_edge = e
                    break

            if existing_edge:
                existing_edge.coverage += 1.0
            else:
                new_edge = Edge(edge_seq)
                new_edge.source = node_left
                new_edge.target = node_right

                node_left.out_edges.append(new_edge)
                node_right.in_edges.append(new_edge)
                self.edges.append(new_edge)

    def compress_graph(self):
        """Сжатие графа: объединение неветвящихся путей."""
        compressed = True
        while compressed:
            compressed = False

            # Находим первый попавшийся узел, который можно сжать
            # Условие: ровно 1 вход, ровно 1 выход, и он не является самозацикленным
            node_to_compress = None
            for node in self.nodes.values():
                if len(node.in_edges) == 1 and len(node.out_edges) == 1:
                    if node.in_edges[0].source != node and node.out_edges[0].target != node:
                        node_to_compress = node
                        break

            if node_to_compress:
                e_in = node_to_compress.in_edges[0]
                e_out = node_to_compress.out_edges[0]

                # Формируем новую последовательность: левая + правая (без учета перекрывающейся части k)
                new_seq = e_in.sequence + e_out.sequence[self.k:]

                # Усредняем покрытие (взвешенное среднее по длине)
                len_in = len(e_in.sequence)
                len_out = len(e_out.sequence)
                new_cov = (e_in.coverage * len_in + e_out.coverage * len_out) / (len_in + len_out)

                new_edge = Edge(new_seq)
                new_edge.coverage = new_cov
                new_edge.source = e_in.source
                new_edge.target = e_out.target

                # Переподключаем source левого ребра
                e_in.source.out_edges.remove(e_in)
                e_in.source.out_edges.append(new_edge)

                # Переподключаем target правого ребра
                e_out.target.in_edges.remove(e_out)
                e_out.target.in_edges.append(new_edge)

                # Удаляем старые связи
                self.edges.remove(e_in)
                self.edges.remove(e_out)
                self.edges.append(new_edge)

                # Удаляем сам узел
                del self.nodes[node_to_compress.kmer]

                compressed = True

        print("Сжатие выполнено.")

    def clean_graph(self, min_coverage: int, max_tip_length: int):
        """Эвристики для очистки графа от ошибок секвенирования."""

        # 1. Удаление ребер с низким покрытием
        edges_to_remove = [e for e in self.edges if e.coverage < min_coverage]
        for e in edges_to_remove:
            e.source.out_edges.remove(e)
            e.target.in_edges.remove(e)
            self.edges.remove(e)

        print(f"Удалено {len(edges_to_remove)} ребер с низким покрытием.")

        # 2. Удаление тупиков (tips)
        # Так как граф уже сжат, тупик — это ребро, ведущее в узел-тупик, и длина этого ребра < max_tip_length
        tip_removed = True
        tips_count = 0
        while tip_removed:
            tip_removed = False
            for e in self.edges:
                # Входящий тупик (в узел e.target ничего больше не входит, и из него ничего не выходит)
                is_out_tip = (len(e.target.out_edges) == 0 and len(e.target.in_edges) == 1)
                # Исходящий тупик
                is_in_tip = (len(e.source.in_edges) == 0 and len(e.source.out_edges) == 1)

                if (is_out_tip or is_in_tip) and len(e.sequence) <= max_tip_length:
                    e.source.out_edges.remove(e)
                    e.target.in_edges.remove(e)
                    self.edges.remove(e)
                    tip_removed = True
                    tips_count += 1
                    break  # Прерываем цикл, чтобы обновить графы безопасно

        print(f"Удалено {tips_count} тупиковых ветвей.")

        # Убираем осиротевшие узлы из словаря
        self.nodes = {k: v for k, v in self.nodes.items() if v.in_edges or v.out_edges}

        # После удаления частей графа могут появиться новые линейные пути, сжимаем еще раз!
        self.compress_graph()
        print("Очистка завершена.")

    def save_to_gfa(self, filepath: str):
        """Сохранение графа в формате GFA v1 для просмотра в Bandage."""
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w') as f:
            f.write("H\tVN:Z:1.0\n")

            edge_id_map = {}
            for i, edge in enumerate(self.edges):
                e_id = f"e{i}"
                edge_id_map[edge] = e_id
                cov_rounded = round(edge.coverage)
                f.write(f"S\t{e_id}\t{edge.sequence}\tKC:i:{cov_rounded}\n")

            for node in self.nodes.values():
                for in_e in node.in_edges:
                    for out_e in node.out_edges:
                        id1 = edge_id_map[in_e]
                        id2 = edge_id_map[out_e]
                        f.write(f"L\t{id1}\t+\t{id2}\t+\t{self.k}M\n")

    def save_contigs_fasta(self, filepath: str):
        with open(filepath, 'w') as f:
            for i, edge in enumerate(self.edges):
                cov = round(edge.coverage, 2)
                f.write(f">contig_{i} length_{len(edge.sequence)} cov_{cov}\n")
                f.write(f"{edge.sequence}\n")


def process_fasta_fastq(filepath: str, format_type: str, k: int) -> DeBruijnGraph:
    graph = DeBruijnGraph(k)
    for record in SeqIO.parse(filepath, format_type):
        graph.add_sequence(str(record.seq))
    return graph


if __name__ == "__main__":
    os.makedirs("outputs", exist_ok=True)

    print("--- Обработка референса ---")
    for k in [15, 31, 55]:
        print(f"Построение графа референса для k={k}")
        graph = process_fasta_fastq("data/ecoli_1k.fna", "fasta", k)
        graph.compress_graph()
        graph.save_to_gfa(f"outputs/ref_k{k}.gfa")

    print("\n--- Обработка прочтений ---")
    k_reads = 31
    print(f"Построение сырого графа чтений для k={k_reads}")
    graph_reads = process_fasta_fastq("data/ecoli_reads.fastq", "fastq", k_reads)
    graph_reads.compress_graph()
    graph_reads.save_to_gfa(f"outputs/reads_raw_k{k_reads}.gfa")

    print("\nОчистка графа чтений")
    # Эвристика: удаляем покрытие < 3 и тупики короче 2 * k
    graph_reads.clean_graph(min_coverage=3, max_tip_length=2 * k_reads)

    graph_reads.save_to_gfa(f"outputs/reads_cleaned_k{k_reads}.gfa")
    graph_reads.save_contigs_fasta("outputs/reads_contigs.fasta")

    print("\nГотово!")
