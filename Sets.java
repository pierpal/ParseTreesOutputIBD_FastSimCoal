import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeMap;

/**
 * Parses Newick trees and outputs IBD segments
 * @author Pier Palamara <pier@cs.columbia.edu>
 */

public class Sets {

    static boolean setOpen = false;
    static StringBuffer sb = new StringBuffer("");
    public static HashMap<String, Integer> id2int = new HashMap<String, Integer>();
    public static HashMap<Integer, String> int2id = new HashMap<Integer, String>();

    public Sets() {
    }

    public static String[][] getSets(String tree) throws IOException {
        /* LOAD THE TREE */
//        System.out.println(tree);
        PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        NexusPhylogeniesParser parser = new NexusPhylogeniesParser();
        Phylogeny[] phylogenies1 = factory.create(new StringBuffer(tree), parser);
        Phylogeny phyl = phylogenies1.clone()[0];
        int N = phyl.getExternalNodes().size();

        String[][] common = new String[N][N];

        TreeMap<Integer, SimpleNode> simpleTree = new TreeMap<Integer, SimpleNode>();
        LinkedList todo = new LinkedList<SimpleNode>();
        Set<PhylogenyNode> external = phyl.getExternalNodes();
        int count=0;
        for (PhylogenyNode ext : external) {
            int ID;
            if (id2int.containsKey(ext.getNodeName())) {
                ID = id2int.get(ext.getNodeName());
            }
            else {
                ID = count++;
//                System.out.println(ID + " is now " + ext.getNodeName());
                id2int.put(ext.getNodeName(), ID);
                int2id.put(ID, ext.getNodeName());
//                System.out.println(ID);
            }
            SimpleNode node = new SimpleNode(ext.getNodeName(), ext.getNodeId(), 0.);
//            System.out.println("ID " + ext.getNodeId() + " " + ext.getDistanceToParent());
            node.descendants.add(ID);
            simpleTree.put(ext.getNodeId(), node);
            todo.push(node);
        }
        // todo is a linked list where I initially push all leaves

        PhylogenyNode root = phyl.getRoot();
        double height = 0.;
        int counter = 0;
        while (todo.size() > 0) { //for each node in the todo list, climb up the tree and build the descendant sets.
            counter++;
            SimpleNode node = (SimpleNode) todo.pollFirst();
            PhylogenyNode node_i = phyl.getNode(node.ID);
            if (node_i.getNodeId() == root.getNodeId()) {
                continue;
            }
//            System.out.println("fetched " + node_i.getNodeId() + " " + node_i.getParent().getNodeId() + " " + node_i.getDistanceToParent());
            PhylogenyNode parent = node_i.getParent();
            Integer parentID = parent.getNodeId();
            double parentHeight = node_i.getDistanceToParent() + node.height;
            String str = String.valueOf(parentHeight); // the name is HEIGHT.ANCESTOR_ID
            String parentNameGen = str;
//            System.out.println("parent of " + node.name + " is " + parentNameGen + " at " + node.height + " yeah " + str + " dtp " + node_i.getDistanceToParent());
            SimpleNode p;
            if (simpleTree.containsKey(parentID)) {
                p = simpleTree.get(parentID);
            } else {
//                System.out.println(parentNameGen + " " + parentID + " " + parentHeight);
                p = new SimpleNode(parentNameGen, parentID, parentHeight);
            }
            p.descendants.addAll(node.descendants);
            p.children.add(node);
            simpleTree.put(parentID, p);
            todo.addLast(p);
//            System.out.println("count " + counter);
        }

        todo.clear();
        counter = 0;
        todo.addLast(simpleTree.get(root.getNodeId()));
        while (todo.size() > 0) {
            counter++;
            SimpleNode node = (SimpleNode) todo.pollFirst();

            for (SimpleNode child : node.children) {
//                System.out.print(child.name + ", ");
                todo.addLast(child);
            }
//            System.out.println("for " + node.name);


            Object[] sets = node.children.toArray();
            for (int i = 0; i < sets.length; i++) {
                for (int j = i + 1; j < sets.length; j++) { // FOR ALL PAIRS OF CHILDREN
                    for (Integer n1 : ((SimpleNode) sets[i]).descendants) {
                        for (Integer n2 : ((SimpleNode) sets[j]).descendants) { // FOR ALL PAIRS OF DESCENDANTS
//                            System.out.println("For " + n1 + " and " + n2 + " parent: " + node.name);
                            if (n2 > n1) {
                                common[n1][n2] = node.name;
                            } else {
                                common[n2][n1] = node.name;
                            }
                        }
                    }
                }
            }
        }
//        for (int i = 0; i < N; i++) {
//            for (int j = 0; j < N; j++) {
//                System.out.println(i + " " + j + " " + common[i][j]);
//            }
//        }
        return common;
    }

    public static void preorder(PhylogenyNode N, double D, double height) {

        if (N.getDistanceToParent() >= 0) {
            D += N.getDistanceToParent();
        }

        int heightFromGround = (int) (height - D);
        ArrayList<PhylogenyNode> descend = N.getDescendants();
        System.out.println(heightFromGround);

        for (int i = 0; i < descend.size(); i++) {
            preorder(descend.get(i), D, height);
        }
    }
}

class SimpleNode {

    String name;
    int ID;
    double height;
    HashSet<Integer> descendants;
    SimpleNode parent;
    HashSet<SimpleNode> children;

    public SimpleNode(String name, int ID, double h) {
        this.name = name;
        this.ID = ID;
        this.height = h;
        this.descendants = new HashSet<Integer>();
        this.children = new HashSet<SimpleNode>();
    }
}
