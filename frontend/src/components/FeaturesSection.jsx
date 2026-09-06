import React, { useEffect, useRef, useState } from 'react';
import { 
  Zap, Filter, Dna, Brain, Search, 
  GitBranch, ScatterChart, Gem, PieChart, X 
} from 'lucide-react';
import './FeaturesSection.css';

const FeaturesSection = ({ onNavigate }) => {
  const sectionRef = useRef(null);
  const [selectedFeature, setSelectedFeature] = useState(null);
  const [isModalOpen, setIsModalOpen] = useState(false);

  const featuresData = [
    {
      id: 'fastp',
      icon: Zap,
      iconImage: '/assets/images/tools/quality-control.png',
      title: 'Quality Control',
      description: 'Lightning-fast preprocessing with fastp for pristine data',
      link: 'tools/fastp',
      // Extended information for modal
      brief: {
        overview: 'Fastp is an ultra-fast all-in-one FASTQ preprocessor that provides quality control, adapter trimming, quality filtering, and per-read quality pruning.',
        keyFeatures: [
          'Automatic adapter detection and removal',
          'Quality filtering and base correction',
          'HTML and JSON reports for quality metrics',
          'Support for both single-end and paired-end data'
        ],
        benefits: [
          'Process millions of reads in minutes',
          'Improve downstream analysis accuracy',
          'Reduce storage requirements by removing low-quality data'
        ]
      }
    },
    {
      id: 'dada2',
      icon: Filter,
      iconImage: '/assets/images/tools/denoising.png',
      title: 'Denoising & ASVs',
      description: 'Gold-standard sequence variant detection with DADA2',
      link: 'tools/dada2',
      brief: {
        overview: 'DADA2 is a state-of-the-art denoising algorithm that infers exact amplicon sequence variants (ASVs) from high-throughput amplicon sequencing data.',
        keyFeatures: [
          'Model-based approach for error correction',
          'Single-nucleotide resolution',
          'Chimera detection and removal',
          'Paired-end read merging'
        ],
        benefits: [
          'Higher resolution than traditional OTU clustering',
          'Reproducible results across studies',
          'Better detection of rare variants'
        ]
      }
    },
    {
      id: 'kraken2',
      icon: Dna,
      iconImage: '/assets/images/tools/taxonomic.png',
      title: 'Taxonomic Classification',
      description: 'Ultra-fast species identification with Kraken2',
      link: 'tools/kraken2',
      brief: {
        overview: 'Kraken2 is a taxonomic classification system that uses exact k-mer matches to achieve high accuracy and fast classification speeds.',
        keyFeatures: [
          'Classify millions of reads per minute',
          'Custom database support',
          'Confidence scoring for each classification',
          'Memory-efficient probabilistic data structures'
        ],
        benefits: [
          'Real-time species identification',
          'High sensitivity for detecting low-abundance organisms',
          'Scalable to massive metagenomic datasets'
        ]
      }
    },
    {
      id: 'dnabert',
      icon: Brain,
      iconImage: '/assets/images/tools/ai-analysis.png',
      title: 'AI-Powered Analysis',
      description: 'Deep learning with DNABERT & SimCLR',
      link: 'tools/ai-analysis',
      brief: {
        overview: 'Leverage state-of-the-art deep learning models DNABERT and SimCLR for advanced sequence analysis and pattern recognition in environmental DNA.',
        keyFeatures: [
          'Pre-trained transformer models for DNA sequences',
          'Self-supervised learning for unlabeled data',
          'Feature extraction for downstream tasks',
          'Transfer learning capabilities'
        ],
        benefits: [
          'Discover hidden patterns in sequence data',
          'Improve classification accuracy',
          'Handle novel or poorly characterized sequences'
        ]
      }
    },
    {
      id: 'faiss',
      icon: Search,
      iconImage: '/assets/images/tools/novelty.png',
      title: 'Novelty Detection',
      description: 'Discover unknown species with FAISS & HDBSCAN',
      link: 'tools/novelty',
      brief: {
        overview: 'Combine Facebook AI Similarity Search (FAISS) with HDBSCAN clustering to identify potentially novel species in your environmental samples.',
        keyFeatures: [
          'Efficient similarity search in high-dimensional space',
          'Density-based clustering for outlier detection',
          'Scalable to billions of sequences',
          'Interactive visualization of clusters'
        ],
        benefits: [
          'Identify candidates for new species discovery',
          'Detect contamination or unusual sequences',
          'Prioritize samples for further investigation'
        ]
      }
    },
    {
      id: 'epang',
      icon: GitBranch,
      iconImage: '/assets/images/tools/phylogenetic.png',
      title: 'Phylogenetic Placement',
      description: 'Evolutionary analysis with EPA-ng',
      link: 'tools/phylogenetics',
      brief: {
        overview: 'EPA-ng (Evolutionary Placement Algorithm) places query sequences into existing phylogenetic trees to understand evolutionary relationships.',
        keyFeatures: [
          'Maximum likelihood placement algorithm',
          'Support for large reference trees',
          'Uncertainty quantification',
          'Parallel processing capabilities'
        ],
        benefits: [
          'Understand evolutionary context of sequences',
          'Identify closest known relatives',
          'Support biodiversity and conservation studies'
        ]
      }
    },
    {
      id: 'umap',
      icon: ScatterChart,
      iconImage: '/assets/images/tools/visualization.png',
      title: 'Data Visualization',
      description: 'Interactive clustering with UMAP',
      link: 'tools/visualization',
      brief: {
        overview: 'UMAP (Uniform Manifold Approximation and Projection) creates stunning visualizations of complex biodiversity data in 2D or 3D space.',
        keyFeatures: [
          'Preserve both local and global data structure',
          'Interactive exploration of clusters',
          'Customizable distance metrics',
          'Integration with diversity indices'
        ],
        benefits: [
          'Identify community patterns at a glance',
          'Detect outliers and anomalies',
          'Communicate findings effectively'
        ]
      }
    },
    {
      id: 'diamond',
      icon: Gem,
      iconImage: '/assets/images/tools/protein.png',
      title: 'Protein Homology',
      description: 'Remote similarity search with DIAMOND',
      link: 'tools/protein',
      brief: {
        overview: 'DIAMOND accelerates BLAST-like protein searches by up to 20,000x, enabling functional annotation of environmental sequences.',
        keyFeatures: [
          'Ultra-fast protein alignment',
          'Sensitive homology detection',
          'Support for multiple databases',
          'Frameshift-aware alignment'
        ],
        benefits: [
          'Predict protein functions from DNA sequences',
          'Identify metabolic pathways in communities',
          'Link taxonomy to function'
        ]
      }
    },
    {
      id: 'diversity',
      icon: PieChart,
      iconImage: '/assets/images/tools/diversity.png',
      title: 'Diversity Analysis',
      description: 'Shannon, Simpson & Bray-Curtis metrics',
      link: 'tools/diversity',
      brief: {
        overview: 'Comprehensive diversity analysis toolkit featuring alpha and beta diversity metrics to quantify and compare biodiversity across samples.',
        keyFeatures: [
          'Alpha diversity: Shannon, Simpson, Chao1',
          'Beta diversity: Bray-Curtis, UniFrac, Jaccard',
          'Rarefaction curves and sample size estimation',
          'Statistical testing for significance'
        ],
        benefits: [
          'Quantify biodiversity changes over time',
          'Compare communities across environments',
          'Support conservation decision-making'
        ]
      }
    }
  ];

  useEffect(() => {
    const observerOptions = {
      threshold: 0.1,
      rootMargin: '0px 0px -50px 0px'
    };

    const observer = new IntersectionObserver((entries) => {
      entries.forEach(entry => {
        if (entry.isIntersecting) {
          entry.target.classList.add('animate-in');
        }
      });
    }, observerOptions);

    const cards = sectionRef.current?.querySelectorAll('.feature-card');
    cards?.forEach(card => observer.observe(card));

    return () => observer.disconnect();
  }, []);

  const handleCardClick = (feature) => {
    setSelectedFeature(feature);
    setIsModalOpen(true);
  };

  const closeModal = () => {
    setIsModalOpen(false);
    setTimeout(() => setSelectedFeature(null), 300); // Clear after animation
  };

  const handleLearnMore = () => {
    if (onNavigate && selectedFeature) {
      onNavigate(selectedFeature.link);
      closeModal();
    }
  };

  return (
    <>
      <section className="features-section" ref={sectionRef}>
        <div className="features-container">
          <div className="features-header">
            <span className="features-overline">COMPREHENSIVE TOOLKIT</span>
            <h2 className="features-main-title">
              Everything you need for modern eDNA discovery
            </h2>
            <p className="features-subtitle">
              Industry-leading tools and AI models built for next-generation environmental genomics
            </p>
          </div>
          
          <div className="features-grid">
            {featuresData.map((feature, index) => {
              const Icon = feature.icon;
              return (
                <div 
                  key={feature.id}
                  className="feature-card" 
                  data-tool={feature.id}
                  style={{ animationDelay: `${index * 0.05}s` }}
                  onClick={() => handleCardClick(feature)}
                >
                  <div className="feature-card-content">
                    <div className="feature-icon-wrapper">
                      {feature.iconImage ? (
                        <img 
                          src={feature.iconImage} 
                          alt={`${feature.title} icon`}
                          className="feature-icon-image"
                        />
                      ) : (
                        <Icon size={24} className="feature-icon" />
                      )}
                    </div>
                    <div className="feature-text">
                      <h3 className="feature-title">{feature.title}</h3>
                      <p className="feature-description">{feature.description}</p>
                    </div>
                  </div>
                </div>
              );
            })}
          </div>

          <div className="features-cta">
            <button 
              className="btn-explore-all"
              onClick={() => onNavigate && onNavigate('tools')}
            >
              Explore all tools
            </button>
          </div>
        </div>
      </section>

      {/* Modal */}
      {selectedFeature && (
        <div className={`feature-modal-overlay ${isModalOpen ? 'open' : ''}`} onClick={closeModal}>
          <div className="feature-modal" onClick={(e) => e.stopPropagation()}>
            <button className="modal-close" onClick={closeModal}>
              <X size={24} />
            </button>
            
            <div className="modal-header">
              <div className="modal-icon-wrapper">
                {selectedFeature.iconImage ? (
                  <img 
                    src={selectedFeature.iconImage} 
                    alt={`${selectedFeature.title} icon`}
                    className="modal-icon-image"
                  />
                ) : (
                  <selectedFeature.icon size={32} className="modal-icon" />
                )}
              </div>
              <h3 className="modal-title">{selectedFeature.title}</h3>
            </div>

            <div className="modal-content">
              <p className="modal-overview">{selectedFeature.brief.overview}</p>
              
              <div className="modal-section">
                <h4 className="modal-section-title">Key Features</h4>
                <ul className="modal-feature-list">
                  {selectedFeature.brief.keyFeatures.map((feature, index) => (
                    <li key={index}>{feature}</li>
                  ))}
                </ul>
              </div>

              <div className="modal-section">
                <h4 className="modal-section-title">Benefits</h4>
                <ul className="modal-benefit-list">
                  {selectedFeature.brief.benefits.map((benefit, index) => (
                    <li key={index}>{benefit}</li>
                  ))}
                </ul>
              </div>
            </div>

            <div className="modal-footer">
              <button className="btn-modal-secondary" onClick={closeModal}>
                Close
              </button>
              <button className="btn-modal-primary" onClick={handleLearnMore}>
                Learn More
              </button>
            </div>
          </div>
        </div>
      )}
    </>
  );
};

export default FeaturesSection;